//
//  opt_bfgs.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::BFGSMin<TT>::init(Vec<TT> vstate)
{
    IM_CHECK_VALID(vstate);
    IM_CHECK_ARGS(vstate.rows()>0);
    
    m_vstate = vstate;
    m_early_exit = false;
    m_hessianok = false;
    m_startup = false;
    m_fx = (TT)0;
    m_delta_fx = (TT)0;
    m_delta_x = (TT)0;
    
    m_iterations = 0;
    
    int const d = dims();
    
    m_linemin.init(this, d, m_params.line_min_eps);
    
    m_vprevstate.resize(d);
    m_vgrad.resize(d);
    m_vprevgrad.resize(d);
    m_vdelta_state.resize(d);
    m_vdelta_grad.resize(d);
    m_mH.resize(d,d);
    m_vHg.resize(d);
    m_vgH.resize(d);
    eval_init();
}

template <typename TT>
bool im::BFGSMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    bool complete = false;
    
    if(m_startup)
    {
        m_mH.set_identity();
        
        m_fx = eval_fx(m_vstate);
        eval_dfx(m_vgrad, m_vstate);

        m_startup = false;
    }
    else
    {
        eval_dfx(m_vgrad, m_vstate);
        
        core_block_sub_elements(m_vdelta_state.view(), m_vstate.view(), m_vprevstate.view());
        core_block_sub_elements(m_vdelta_grad.view(), m_vgrad.view(), m_vprevgrad.view());
        
        TT dx_dot_dg = m_vdelta_state.dot_product(m_vdelta_grad);
        
        if(m_hessianok==false)
        {
            TT dg_dot_dg = m_vdelta_grad.magnitude_squared();
            if(dg_dot_dg > TypeProperties<TT>::epsilon())
            {
                TT t = dx_dot_dg / dg_dot_dg;
                if(t<(TT)0.01)
                    t = (TT)0.01;
                else if(t>(TT)100.0)
                    t = (TT)100.0;
                m_mH.diag() = t;
                m_hessianok = true;
            }
        }
        
        core_block_blas_gemv(m_vHg.view(), m_mH.view(), m_vdelta_grad.view(), (TT)1, (TT)0, TransMode_N);
        core_block_blas_gemv(m_vgH.view(), m_mH.view(), m_vdelta_grad.view(), (TT)1, (TT)0, TransMode_T);
        TT gHg = m_vdelta_grad.dot_product(m_vHg);
        
        if(gHg < std::numeric_limits<TT>::infinity() && dx_dot_dg < std::numeric_limits<TT>::infinity() && dx_dot_dg != (TT)0)
        {
            // update H
            core_block_blas_ger(m_mH.view(), m_vdelta_state.view(), m_vdelta_state.view(), ((TT)1 + gHg / dx_dot_dg)/dx_dot_dg, TransMode_N);
            core_block_blas_ger(m_mH.view(), m_vdelta_state.view(), m_vgH.view(), -(TT)1/dx_dot_dg, TransMode_N);
            core_block_blas_ger(m_mH.view(), m_vHg.view(), m_vdelta_state.view(), -(TT)1/dx_dot_dg, TransMode_N);
        }
        else
        {
            m_mH.set_identity();
            m_hessianok = false;
        }
    }

    m_vprevstate.copy_from(m_vstate);
    m_vprevgrad.copy_from(m_vgrad);
    core_block_blas_gemv(m_vdelta_grad.view(), m_mH.view(), m_vgrad.view(), (TT)-1, (TT)0, TransMode_N);
    
    TT fxstart = m_fx;
    m_linemin.linemin(m_vstate, m_vdelta_grad, m_params.bracket_max);
    m_fx = m_linemin.fxmin();
    m_delta_x = m_linemin.xmin() * m_vdelta_grad.magnitude_squared();
    m_delta_fx = m_fx - fxstart;
    
    if((TT)2*(std::abs(m_delta_fx)) <= m_params.termination_ratio * (std::abs(m_fx) + std::abs(fxstart) + TypeProperties<TT>::epsilon()))
        complete = true;
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    m_iterations++;
    if(m_iterations > m_params.iterations_max)
        complete = true;
    
    return m_early_exit || complete;
}

template class im::BFGSMin<float>;
template class im::BFGSMin<double>;
