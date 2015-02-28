//
//  decomp_tridiag.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::MatrixDecompTridiag<TT>::compute(MtxView<TT> const &mavA)
{
    IM_CHECK_VALID(mavA);
    IM_CHECK_MATRIX_SQUARE(mavA);

    int size = mavA.rows();
    IM_CHECK_ARGS(size>0);
    
    m_hh.resize(size, size);
    m_vtau.resize(size-1);

    core_block_copy_lower_tri(m_hh.view(), mavA, true);
    core_block_copy_upper_tri(m_hh.view(), m_hh.view().t(), false); // save A
    
    compute_in_place(m_hh.view(), m_vtau.view());
}

template <typename TT>
void im::MatrixDecompTridiag<TT>::compute_in_place(MtxView<TT> mavA, VecView<TT> vTau)
{
    IM_CHECK_VALID(mavA);
    IM_CHECK_MATRIX_SQUARE(mavA);
    IM_CHECK_VALID(vTau);
    
    int size = mavA.rows();
    IM_CHECK_ARGS(size>0);
    IM_CHECK_ARGS(vTau.rows()==size-1);
    
    if(size<3)
    {
        // Matrix is already in tridiag form, so Q = I, T = A.
        vTau = (TT)0;
        return;
    }
    
    int i;
    for(i=0; i<size-2; i++)
    {
        int const block_size = size - (i+1); // number of elements below the diagonal
        
        VecView<TT> v = mavA.col(i).block(i+1, block_size); // column vector below the diagonal v
        TT tau = core_compute_householder_vector(v);
        
        // Apply transformation H' A H to remaining columns
        if(tau!=(TT)0)
        {
            MtxView<TT> m = mavA.block(i+1, i+1, block_size, block_size);
            VecView<TT> x = vTau.block(i, block_size);
            
            TT g = v(0);
            v(0) = (TT)1;
            
            // x = tau A v
            core_block_blas_symv(x, m, v, tau, (TT)0, TriMode_L);
            
            // x = x - (1/2) tau (x' v) v
            TT alpha = -(TT)0.5 * tau * core_block_reduce_multiply_add(x, v);
            core_block_add_scaled(x, x, v, alpha);
            
            // A = A - v w' - w v'
            for(int k=0; k<block_size; k++)
            {
                TT t1 = -x(k);
                TT t2 = -v(k);
                for(int j=0; j<=k; j++)
                    m(k,j) += t1 * v(j) + t2 * x(j);
            }
            
            v(0) = g;
        }
        
        vTau(i) = tau;
    }
}

// Gets the matrix Q
template <typename TT>
im::Mtx<TT> im::MatrixDecompTridiag<TT>::matrix_Q() const
{
    int size = m_hh.rows();
    Mtx<TT> mQ(size,size);
    mQ.set_identity();
    
    if(size>=3)
    {
        for(int i=size-2; i-- > 0; )
        {
            int const block_size = size - (i+1);
            core_apply_householder_vector(mQ.view().block(i+1, i+1, block_size, block_size),
                                          m_vtau(i),
                                          m_hh.view().col(i).block(i+1, block_size));
        }
    }
    
    return mQ;
}

// Gets the matrix T
template <typename TT>
im::Mtx<TT> im::MatrixDecompTridiag<TT>::matrix_T() const
{
    int size = m_hh.rows();
    Mtx<TT> mT(size,size);
    mT = (TT)0;
    
    mT.diag().copy_from(m_hh.diag());
    if(size>1)
    {
        mT.diag(-1).copy_from(m_hh.diag(-1));
        mT.diag(1).copy_from(m_hh.diag(-1));
    }
    
    return mT;
}

template class im::MatrixDecompTridiag<float>;
template class im::MatrixDecompTridiag<double>;

