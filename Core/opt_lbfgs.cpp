//
//  opt_lbfgs.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/17/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::LBFGSMin<TT>::init(Vec<TT> vstate)
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
    
    int const d = dims();
    
    m_linemin.init(this, d, m_params.line_min_eps, true);
    
    m_vprevstate.resize(d);
    m_vgrad.resize(d);
    m_vprevgrad.resize(d);
    m_vprevdir.resize(d);
    
    m_history_length = m_params.history_length;
    if(m_history_length<1)
        m_history_length = d;
    
    eval_init();
}

template <typename TT>
bool im::LBFGSMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    bool complete = false;
    
    if(m_startup)
    {
        m_fx = eval_fx(m_vstate);
        eval_dfx(m_vgrad, m_vstate);
        
        // search downhill
        core_block_neg(m_vprevdir.view(), m_vgrad.view());
        
        m_startup = false;
    }
    else
    {
        eval_dfx(m_vgrad, m_vstate);
        
        // search downhill
        core_block_neg(m_vprevdir.view(), m_vgrad.view());
        
        LBFGSItem<TT> item;
        item.vs = m_vstate - m_vprevstate;
        item.vy = m_vgrad - m_vprevgrad;
        TT s_dot_y = item.vs.dot_product(item.vy);
        
        if(std::abs(s_dot_y) > TypeProperties<TT>::epsilon())
        {
            item.rho = (TT)1 / s_dot_y;
            m_queue.push_back(item); // add new item
        }
        else
            m_queue.clear(); // discard everything and start again
        
        if(m_queue.size()>0)
        {
            // compute better direction than down the gradient
            m_alpha.resize(m_queue.size());
            for(int i = (int)m_queue.size() - 1; i>=0; i--)
            {
                m_alpha[i] = m_queue[i].rho * m_queue[i].vs.dot_product(m_vprevdir);
                core_block_blas_axpy(m_vprevdir.view(), m_queue[i].vy.view(), -m_alpha[i]); // pd -= alpha*vy
            }

            TT h0 = (TT)1 / (m_queue.back().rho * m_queue.back().vy.magnitude_squared());
            h0 = core_range_limit(h0, (TT)0.001, (TT)1000.0);
            core_block_scale(m_vprevdir.view(), m_vprevdir.view(), h0);
            
            for(int i=0; i<(int)m_queue.size(); i++)
            {
                TT beta = m_queue[i].rho * m_queue[i].vy.dot_product(m_vprevdir);
                core_block_blas_axpy(m_vprevdir.view(), m_queue[i].vs.view(), m_alpha[i] - beta);
            }
        }
    }

    // throw away oldest
    if(m_queue.size() > m_history_length)
        m_queue.pop_front();
    
    // update
    m_vprevstate.copy_from(m_vstate);
    m_vprevgrad.copy_from(m_vgrad);
    
    TT fxstart = m_fx;
    TT dirmag = m_vprevdir.magnitude();
    if(dirmag > TypeProperties<TT>::epsilon())
    {
        m_linemin.linemin(m_vstate, m_vprevdir, m_params.bracket_max);
        m_fx = m_linemin.fxmin();
        m_delta_x = std::abs(m_linemin.xmin()) * dirmag;
        m_delta_fx = m_fx - fxstart;
        
        if((TT)2*(std::abs(m_delta_fx)) <= m_params.termination_ratio * (std::abs(m_fx) + std::abs(fxstart) + TypeProperties<TT>::epsilon()))
            complete = true;
    }
    else
    {
        m_delta_x = (TT)0;
        m_delta_fx = (TT)0;
        complete = true;
    }
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    m_iterations++;
    if(m_iterations > m_params.iterations_max)
        complete = true;
    
    return m_early_exit || complete;
}

template class im::LBFGSMin<float>;
template class im::LBFGSMin<double>;
