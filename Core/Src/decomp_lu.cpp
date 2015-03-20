//
//  decomp_lu.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
void im::MatrixDecompLU<TT>::compute(MtxView<TT> const &mavA)
{
    IM_CHECK_VALID(mavA);
    IM_CHECK_MATRIX_SQUARE(mavA);
    
    int size = mavA.rows();
    
    m_matLU.resize(size, size);
    m_matLU.copy_from(mavA);
    
    m_perm.resize(size);
    
    for(int j=0; j<size-1; j++)
    {
        MtxLoc loc = m_matLU.block(j,j,size - j,1).location_of_max_abs();
        loc.row += j; // I went crazy for ages because I missed this line out!
        
        if(loc.row != j)
        {
            core_block_exchange(m_matLU.view().row(j), m_matLU.view().row(loc.row));
            m_perm.swap(j, loc.row);
        }
        
        TT ajj = m_matLU(j,j);
        
        if(ajj != (TT)0)
        {
            for(int i = j+1; i<size; i++)
            {
                TT aij = m_matLU(i,j) / ajj;
                m_matLU(i,j) = aij;
                
                for(int k = j+1; k<size; k++)
                    m_matLU(i,k) -= aij * m_matLU(j,k);
            }
        }
    }
}

// Solve Ax = y
template <typename TT>
im::Mtx<TT> im::MatrixDecompLU<TT>::solve(MtxView<TT> const &mavy) const
{
    IM_CHECK_VALID(mavy);
    IM_CHECK_ARGS(mavy.rows()==m_matLU.rows());
    
    int size = mavy.rows();
    int numcols = mavy.cols();
    
    Mtx<TT> mtxrtn(size,numcols);
    
    // apply permutation to y: Py = LUx
    for(int i=0; i<size; i++)
        mtxrtn.view().row(i).copy_from(mavy.row(m_perm(i)));

    // solve Pb = Lc
    core_block_blas_trsm(mtxrtn.view(), m_matLU.view(), (TT)1, SideMode_L, TriMode_L, TransMode_N, DiagMode_U);

    // solve c = Ux
    core_block_blas_trsm(mtxrtn.view(), m_matLU.view(), (TT)1, SideMode_L, TriMode_U, TransMode_N, DiagMode_N);
    
    return mtxrtn;
}

template <typename TT>
im::Vec<TT> im::MatrixDecompLU<TT>::solve(VecView<TT> const &vvy) const
{
    IM_CHECK_VALID(vvy);
    IM_CHECK_ARGS(vvy.rows()==m_matLU.rows());
    
    int size = vvy.rows();
    
    Vec<TT> vrtn(size);
    
    // apply permutation to y: Py = LUx
    for(int i=0; i<size; i++)
        vrtn(i) = vvy(m_perm(i));
    
    // solve Pb = Lc
    core_block_blas_trsv(vrtn.view(), m_matLU.view(), TriMode_L, TransMode_N, DiagMode_U);
    
    // solve c = Ux
    core_block_blas_trsv(vrtn.view(), m_matLU.view(), TriMode_U, TransMode_N, DiagMode_N);
    
    return vrtn;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompLU<TT>::inverse() const
{
    int size = m_matLU.rows();
    
    Mtx<TT> matI(size, size);
    matI.set_identity();
    
    return solve(matI.view());
}

template class im::MatrixDecompLU<float>;
template class im::MatrixDecompLU<double>;
template class im::MatrixDecompLU<im::Cf>;
template class im::MatrixDecompLU<im::Cd>;

