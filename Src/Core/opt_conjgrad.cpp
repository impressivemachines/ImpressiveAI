//
//  opt_conjgrad.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::ConjGradientMin<TT>::init(Vec<TT> vstate)
{
    IM_CHECK_VALID(vstate);
    IM_CHECK_ARGS(vstate.rows()>0);
    
    m_vstate = vstate;
    m_early_exit = false;
    m_fx = (TT)0;
    m_delta_fx = (TT)0;
    m_delta_x = (TT)0;
    
    m_iterations = 0;
    
    int const d = dims();
    m_vdir_x.resize(d);
    m_vdir_g.resize(d);
    m_vdir_h.resize(d);
    
    m_linemin.init(this, d, m_params.line_min_eps, true);
    
    m_startup = true;
    
    eval_init();
}

template <typename TT>
bool im::ConjGradientMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    bool complete = false;
    
    if(m_startup)
    {
        m_fx = eval_fx(m_vstate);
        eval_dfx(m_vdir_x, m_vstate);
        
        core_block_neg(m_vdir_x.view(), m_vdir_x.view());
        m_vdir_g.copy_from(m_vdir_x);
        m_vdir_h.copy_from(m_vdir_x);
        
        m_startup = false;
    }
    
    TT fxold = m_fx;
    
    m_linemin.linemin(m_vstate, m_vdir_x, m_params.bracket_max);
    TT fxnew = m_linemin.fxmin();
    m_delta_x = std::abs(m_linemin.xmin()) * m_vdir_x.magnitude();
    
    if((TT)2*(std::abs(fxnew - m_fx)) <= m_params.termination_ratio * (std::abs(fxnew) + std::abs(m_fx) + TypeProperties<TT>::epsilon()))
        complete = true;
    else
    {
        m_fx = fxnew;
        eval_dfx(m_vdir_x, m_vstate);
        
        // test for gradient relatively close to zero
        TT scale = (TT)1 / std::max((TT)1, std::abs(m_fx));
        TT tmax = (TT)0;
        for(int i=0; i<dims(); i++)
        {
            TT t = scale * std::abs(m_vdir_x(i)) * std::max((TT)1, std::abs(m_vstate(i)));
            if(t>tmax)
                tmax = t;
        }
        
        if(tmax < m_params.gradient_ratio)
            complete = true;
        else
        {
            // compute gamma
            TT gamma_denom = m_vdir_g.magnitude_squared();
            if(gamma_denom==(TT)0)
                complete = true;
            else
            {
                TT gamma = m_vdir_x.magnitude_squared();
                if(m_params.update_mode==ConjGradientUpdateModePolakRibiere)
                    gamma += m_vdir_g.dot_product(m_vdir_x);
                    
                gamma /= gamma_denom;
                
                core_block_neg(m_vdir_g.view(), m_vdir_x.view()); // g = -x
                m_vdir_x.copy_from(m_vdir_g); // x = g
                core_block_blas_axpy(m_vdir_x.view(), m_vdir_h.view(), gamma); // x += gamma * h
                m_vdir_h.copy_from(m_vdir_x); // h = x
            }
        }
    }
    
    m_delta_fx = m_fx - fxold;
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    m_iterations++;
    if(m_iterations > m_params.iterations_max)
        complete = true;
    
    return m_early_exit || complete;
}

template class im::ConjGradientMin<float>;
template class im::ConjGradientMin<double>;
