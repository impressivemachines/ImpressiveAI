//
//  permute.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::Permute::matrix_PX_in_place(MtxView<TT> mavX) const
{
    int size = (int)m_permute_map.size();
    IM_CHECK_VALID(mavX);
    IM_CHECK_ARGS(mavX.rows()==size);
    
    for(int i=0; i<size; i++)
    {
        int k = m_permute_map[i];
        
        while(k>i)
            k = m_permute_map[k];
        
        if(k<i)
            continue;
        
        int pk = m_permute_map[k];
        if(pk==i)
            continue;
        
        int ksave = k;
        int pksave = pk;
        
        for(int col = 0; col < mavX.cols(); col++)
        {
            k = ksave;
            pk = pksave;
            
            TT tmp = mavX(i,col);
            while(pk!=i)
            {
                mavX(k,col) = mavX(pk,col);
                k = pk;
                pk = m_permute_map[k];
            }
            mavX(k,col) = tmp;
        }
    }
}

#define INST(TT) template void im::Permute::matrix_PX_in_place(MtxView<TT> mavX) const
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT>
void im::Permute::matrix_PTX_in_place(MtxView<TT> mavX) const
{
    int size = (int)m_permute_map.size();
    IM_CHECK_VALID(mavX);
    IM_CHECK_ARGS(mavX.rows()==size);
    
    for(int i=0; i<size; i++)
    {
        int k = m_permute_map[i];
        
        while(k>i)
            k = m_permute_map[k];
        
        if(k<i)
            continue;
        
        int pk = m_permute_map[k];
        if(pk==i)
            continue;
        
        int ksave = k;
        int pksave = pk;
        
        for(int col = 0; col < mavX.cols(); col++)
        {
            k = ksave;
            pk = pksave;
            
            TT tmp = mavX(k,col);
            while(pk!=i)
            {
                std::swap(tmp, mavX(pk,col));
                k = pk;
                pk = m_permute_map[k];
            }
            mavX(pk,col) = tmp;
        }
    }
}

#define INST(TT) template void im::Permute::matrix_PTX_in_place(MtxView<TT> mavX) const
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

