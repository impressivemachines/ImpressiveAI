//
//  opt_vstepgrad.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/18/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"


template <typename TT>
void im::VariableStepGradientMin<TT>::init(Vec<TT> vstate)
{
    IM_CHECK_VALID(vstate);
    IM_CHECK_ARGS(vstate.rows()>0);
    
    m_vstate = vstate;
    
    m_early_exit = false;
    m_startup = true;
    m_fx = (TT)0;
    m_delta_fx = (TT)0;
    m_delta_x = (TT)0;
    m_iterations = 0;

    int const d = vstate.rows();
    
    m_vtrialstate.resize(d);
    m_vtrialgrad.resize(d);
    m_vgrad.resize(d);
    
    eval_init();
}

template <typename TT>
bool im::VariableStepGradientMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    bool complete = false;

    TT fxstart = m_fx;
    
    if(m_startup)
    {
        m_fx = eval_fx(m_vstate);
        eval_dfx(m_vgrad, m_vstate);
        
        m_startup = false;
    }
    
    TT gradmag = m_vgrad.magnitude();
    
    while(true)
    {
        // try to find a good step forward
        // x' = -g * x
        core_block_blas_axpy(m_vtrialstate.view(), m_vgrad.view(), -m_stepsize / gradmag);
        
        // evaluate the gradient at the new location
        // gtrial = df(x')
        eval_dfx(m_vtrialgrad, m_vtrialstate);
        
        TT trialgradmag = m_vtrialgrad.magnitude();
        
        IM_CHECK_ARGS(false); // under development
        
        // turning factor
        TT dot = m_vtrialgrad.dot_product(m_vgrad) / (gradmag * trialgradmag);
        
        if(dot > m_params.turning_threshold)
        {
            // accept
            m_vstate.copy_from(m_vtrialstate);
            
            // g' = g * momentum + gtrial * (1-momentum)
            core_block_scale(m_vgrad.view(), m_vgrad.view(), m_params.momentum);
            core_block_blas_axpy(m_vgrad.view(), m_vtrialgrad.view(), 1-m_params.momentum);
            
            m_stepsize *= m_params.stepsize_change_on_accept;
            break;
        }
        else
        {
            // reject
            m_stepsize *= m_params.stepsize_change_on_reject;
        }
    }
    
    
    
    
    if((TT)2*(std::abs(m_delta_fx)) <= m_params.termination_ratio * (std::abs(m_fx) + std::abs(fxstart) + TypeProperties<TT>::epsilon()))
        complete = true;
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    m_iterations++;
    if(m_iterations > m_params.iterations_max)
        complete = true;
    
    return m_early_exit || complete;
}

template class im::StochasticMin<float>;
template class im::StochasticMin<double>;
