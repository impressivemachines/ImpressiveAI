//
//  opt_stochastic.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"


template <typename TT>
void im::StochasticMin<TT>::init(Vec<TT> vstate)
{
    IM_CHECK_VALID(vstate);
    IM_CHECK_ARGS(vstate.rows()>0);
    
    m_vstate = vstate;
    m_vgradient.resize(dims());
    
    m_early_exit = false;
    m_startup = true;
    m_fx = (TT)0;
    m_delta_fx = (TT)0;
    m_delta_x = (TT)0;
    m_iterations = 0;
    m_rate = (TT)m_params.rate_init;
    m_momentum = (TT)m_params.momentum_init;
    
    if(m_params.update_mode!=StochasticUpdateModeSimple)
    {
        m_vintegrator.resize(dims());
        m_vintegrator = (TT)0;
    }
    
    eval_init();
}

template <typename TT>
bool im::StochasticMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    bool complete = false;
    
    if(m_iterations % m_params.fx_eval_interval==0)
    {
        TT fx = eval_fx(m_vstate);
        m_delta_fx = (fx - m_fx) / m_params.fx_eval_interval;
        m_fx = fx;
        if(!m_startup && std::abs(m_delta_fx) < m_params.delta_fx_min)
            complete = true;
    }
    
    // update state
    switch(m_params.update_mode)
    {
        StochasticUpdateModeSimple:
        default:
            // x' = x - rate * gradient(x)
            // compute the gradient at x
            eval_dfx(m_vgradient, m_vstate);
            core_block_blas_axpy(m_vstate.view(), m_vgradient.view(), -m_rate);
            if(m_params.compute_delta_x)
                m_delta_x = m_rate * m_vgradient.magnitude();
            break;
            
        StochasticUpdateModeMomentum:
            // compute the gradient at x
            eval_dfx(m_vgradient, m_vstate);
            
            if(m_startup)
            {
                // v' = -rate * gradient
                core_block_scale(m_vintegrator.view(), m_vgradient.view(), -m_rate);
            }
            else
            {
                // v' = momentum * v - rate * gradient(x)
                core_block_scale(m_vintegrator.view(), m_vintegrator.view(), m_momentum); // v *= momentum
                TT rate = m_rate;
                if(m_params.blend_mode)
                    rate *= ((TT)1-m_momentum);
                core_block_blas_axpy(m_vintegrator.view(), m_vgradient.view(), -rate); // v -= rate*gradient
            }
            
            // x' = x + v'
            core_block_add_elements(m_vstate.view(), m_vstate.view(), m_vintegrator.view()); // x += v
            
            if(m_params.compute_delta_x)
                m_delta_x = m_vintegrator.magnitude();
            break;
            
        StochasticUpdateModeAccelerated:
            if(m_startup)
            {
                // compute the gradient at x
                eval_dfx(m_vgradient, m_vstate);
                
                // v' = -rate * gradient
                core_block_scale(m_vintegrator.view(), m_vgradient.view(), -m_rate);
            }
            else
            {
                // v' = momentum * v - rate * gradient(x + momentum * v)
                core_block_scale(m_vintegrator.view(), m_vintegrator.view(), m_momentum); // v *= momentum
                core_block_add_elements(m_vstate.view(), m_vstate.view(), m_vintegrator.view()); // x += v
                
                // compute the gradient at x + momentum * v
                eval_dfx(m_vgradient, m_vstate);
                
                core_block_sub_elements(m_vstate.view(), m_vstate.view(), m_vintegrator.view()); // x -= v
                TT rate = m_rate;
                if(m_params.blend_mode)
                    rate *= ((TT)1-m_momentum);
                core_block_blas_axpy(m_vintegrator.view(), m_vgradient.view(), -rate); // v -= rate * gradient
            }
            
            // x' = x + v'
            core_block_add_elements(m_vstate.view(), m_vstate.view(), m_vintegrator.view()); // x += v
            
            if(m_params.compute_delta_x)
                m_delta_x = m_vintegrator.magnitude();
            break;
    }
    
    // update learning parameters
    if(m_params.update_rate)
    {
        if(m_iterations>0 && (m_iterations % m_params.rate_interval)==0)
            m_rate *= m_params.rate_multiplier;
    }
    
    if(m_params.update_momentum)
        m_momentum = (TT)std::min(1.0 - 0.5/(m_iterations/m_params.momentum_rate + 1.0), m_params.momentum_max);
    
    
    if(m_startup)
    {
        m_delta_fx = (TT)0;
        m_startup = false;
    }
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    m_iterations++;
    if(m_iterations > m_params.iterations_max)
        complete = true;
    
    return m_early_exit || complete;
}

template class im::StochasticMin<float>;
template class im::StochasticMin<double>;
