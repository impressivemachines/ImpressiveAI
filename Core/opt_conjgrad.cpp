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
    
    eval_init();
}

template <typename TT>
bool im::ConjGradientMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    int const d = dims();
    
    bool complete = false;
    
    
    
    
    
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    return m_early_exit || complete;
}

template class im::ConjGradientMin<float>;
template class im::ConjGradientMin<double>;
