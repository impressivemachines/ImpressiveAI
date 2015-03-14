//
//  opt_powell.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::PowellMin<TT>::init(Vec<TT> vstate)
{
    IM_CHECK_VALID(vstate);
    IM_CHECK_ARGS(vstate.rows()>0);
    
    m_vstate = vstate;
    m_early_exit = false;
    m_startup = true;
    m_fx = (TT)0;
    m_delta_fx = (TT)0;
    m_delta_x = (TT)0;
    
    int const d = dims();
    
    m_linemin.init(this, d, m_params.line_min_eps);
    
    m_mdirset.resize(d, d);
    m_mdirset.set_identity();
    m_vdir.resize(d);
    m_vstatesave.resize(d);
    
    m_iterations = 0;
    
    eval_init();
}

template <typename TT>
bool im::PowellMin<TT>::step()
{
    if(m_early_exit)
        return true;
    
    eval_start_step();
    
    int const d = dims();
    
    bool complete = false;
    if(m_startup)
    {
        m_vstatesave.copy_from(m_vstate);
        m_fx = eval_fx(m_vstate);
        m_startup = false;
    }
    
    TT fxstart = m_fx;
    TT maxchange = (TT)0;
    int maxi = 0;
    
    for(int i=0; i<d; i++)
    {
        // line min, moving vstate in the process
        m_linemin.linemin(m_vstate, m_mdirset.col(i), m_params.bracket_max);
        TT fxnew = m_linemin.fxmin();
        TT change = m_fx - fxnew;
        m_fx = fxnew;
        if(change > maxchange)
        {
            maxi = i;
            maxchange = change;
        }
    }
    
    if(fxstart - m_fx <= 0.5f*(std::abs(fxstart) + std::abs(m_fx)) * (TT)m_params.termination_ratio + TypeProperties<TT>::epsilon())
    {
        m_delta_x = (TT)0;
        complete = true;
    }
    else
    {
        // new mean direction = Pn - P0
        core_block_sub_elements(m_vdir.view(), m_vstate.view(), m_vstatesave.view());
    
        m_delta_x = m_vstate.distance(m_vstatesave);
        
        // save Pn as P0 for next step call
        m_vstatesave.copy_from(m_vstate);
        
        // extrapolated point = Pn + (Pn - P0)
        core_block_add_elements(m_vstate.view(), m_vstate.view(), m_vdir.view());
        
        // eval at extrapolated point
        TT fxe = eval_fx(m_vstate);
        
        // restore state
        m_vstate.copy_from(m_vstatesave);
        
        if(fxe < fxstart)
        {
            TT a = fxstart - m_fx - maxchange; // reduction without the largest one
            TT b = fxstart - fxe; // reduction with extrapolation
            TT c = fxstart - (TT)2 * m_fx + fxe;
            if((TT)2 * c * a * a < maxchange * b * b)
            {
                // worth incorpotating the new mean direction into the direction set
                m_linemin.linemin(m_vstate, m_vdir, m_params.bracket_max);
                m_fx = m_linemin.fxmin();
                
                // get rid of the direction with largest reduction
                if(maxi != d-1)
                    m_mdirset.col(maxi).copy_from(m_mdirset.col(d-1));
                
                // replace by scaled mean direction
                core_block_scale(m_mdirset.col(d-1).view(), m_vdir.view(), m_linemin.xmin());
            }
        }
    }
    
    m_delta_fx = m_fx - fxstart;
    
    m_early_exit = eval_end_step() || m_early_exit;
    
    m_iterations++;
    if(m_iterations > m_params.iterations_max)
        complete = true;
    
    return m_early_exit || complete;
}

template class im::PowellMin<float>;
template class im::PowellMin<double>;
