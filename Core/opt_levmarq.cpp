//
//  opt_levmarq.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::LevenbergMarquardt<TT>::init(Vec<TT> vstate)
{
    IM_CHECK_VALID(vstate);
    IM_CHECK_ARGS(vstate.rows()>0);
    
    int const d = vstate.rows();
    m_vstate = vstate;
    m_vresidual.resize(d);
    m_vnewresidual.resize(d);
    m_mJ.resize(d,d);
    m_mJtJ.resize(d,d);
    m_vJte.resize(d);
    m_vJtJdiag.resize(d);
    
    m_early_exit = false;
    m_error = (TT)0;
    m_delta_error = (TT)0;
    m_delta_x = (TT)0;
    
    m_lambda = (TT)m_params.lambda_start;
    
    eval_init();
}

template <typename TT>
bool im::LevenbergMarquardt<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    eval_residual(m_vresidual, m_vstate);
    m_error = m_vresidual.magnitude_squared();
    
    eval_jacobian(m_mJ, m_vstate);
    
    int const d = dims();
    
    // Compute Jte
    core_block_blas_gemv(m_vJte.view(), m_mJ.view(), m_vresidual.view(), (TT)1, (TT)0, TransMode_T);
    
    // Compute JtJ
    for(int i=0; i<d; i++)
        for(int j=0; j<=i; j++)
        {
            TT sum = (TT)0;
            for(int k=0; k<d; k++)
                sum += m_mJ(k,i) * m_mJ(k,j);
            m_mJtJ(i,j) = sum;
            m_mJtJ(j,i) = sum;
        }
    
    // save the diagonal
    for(int i=0; i<d; i++)
        m_vJtJdiag(i) = m_mJtJ(i,i);

    bool complete = false;
    
    // Loop where we try changing lambda to get an improvement
    while(true)
    {
        // augment the normal equations
        if(m_params.use_marquardt_mode)
            for(int i=0; i<d; i++)
                m_mJtJ(i, i) = m_vJtJdiag(i) * ((TT)1 + m_lambda);
        else
            for(int i=0; i<d; i++)
                m_mJtJ(i, i) = m_vJtJdiag(i) + m_lambda;
        
        // Solve J'J delta = J'e to get delta
        m_ldlt.compute(m_mJtJ.view());
       
        Vec<TT> vnewstate = m_ldlt.solve(m_vJte.view());
        
        if(m_params.update_fraction < 1)
            core_block_scale(vnewstate.view(), vnewstate.view(), (TT)m_params.update_fraction);
        
        core_block_sub_elements(vnewstate.view(), m_vstate.view(), vnewstate.view());
        
        eval_residual(m_vnewresidual, vnewstate);
        TT newerror = m_vnewresidual.magnitude_squared();
        
        bool better;
        if(newerror < m_error)
        {
            // better so decrease lambda
            m_lambda /= m_params.lambda_ratio;
            better = true;
        }
        else
        {
            m_lambda *= m_params.lambda_ratio;
            better = false;
        }
        
        // termination check
        TT error_ratio_sq = core_block_reduce_squared_distance(m_vresidual.view(), m_vnewresidual.view()) / (m_error + TypeProperties<TT>::epsilon());
        if(error_ratio_sq < m_params.error_ratio_min * m_params.error_ratio_min)
            complete = true;

        if(m_lambda > m_params.lambda_max)
        {
            if(m_params.fail_on_lambda_max)
            {
                // failure case of getting nowhere
                m_early_exit = true;
            }
            else
                m_lambda = m_params.lambda_max;
            
            complete = true;
        }
        
        if(better)
        {
            m_delta_x = (vnewstate - m_vstate).magnitude();
            m_delta_error = newerror - m_error;
            m_error = newerror;
            m_vstate.copy_from(vnewstate);
        }
        
        if(better || complete)
            break;
    }

    m_early_exit = eval_end_step() || m_early_exit;
    
    return m_early_exit || complete;
}

template <typename TT>
void im::LevenbergMarquardt<TT>::eval_residual_diff(Vec<TT> &vdiff, Vec<TT> &vx, int jcol, TT delta, Vec<TT> const &vfx)
{
    vx(jcol) += delta;
    
    eval_residual(vdiff, vx);
    vdiff -= vfx;
}

template <typename TT>
void im::LevenbergMarquardt<TT>::eval_jacobian(Mtx<TT> &mjacobian, Vec<TT> &vx)
{
    int d = dims();
    
    // residual is already pre-computed for vx in m_vresidual, so no need to repeat here
    // loop over columns of J
    for(int j=0; j<d; j++)
    {
        // find increment
        TT dx = std::abs(vx(j) * m_params.deriv_delta_ratio);
        if(dx < m_params.deriv_delta_min)
            dx = m_params.deriv_delta_min;
        
        // get a column of the jacobian
        TT oldx = vx(j);
        eval_residual_diff(m_vnewresidual, vx, j, dx, m_vresidual);
        vx(j) = oldx;
        
        // update column
        for(int i=0; i<d; i++)
            mjacobian(i,j) = m_vnewresidual(i) / dx;
    }
}

template class im::LevenbergMarquardt<float>;
template class im::LevenbergMarquardt<double>;

