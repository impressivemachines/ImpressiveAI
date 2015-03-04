//
//  decomp_cholesky.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::MatrixDecompLDLT<TT>::compute(MtxView<TT> const &mavA)
{
    IM_CHECK_VALID(mavA);
    IM_CHECK_MATRIX_SQUARE(mavA);
    
    m_matLDLT.resize(mavA.rows(), mavA.cols());
    
    core_block_copy_lower_tri(m_matLDLT.view(), mavA, true);
    
    int size = mavA.rows();
    
    int i, k, row;
    
    for(row = 0; row<size; row++)
    {
        for(i=0; i<=row-1; i++)
            m_matLDLT(i,row) = m_matLDLT(row,i) * m_matLDLT(i,i);
        
        TT tmp = m_matLDLT(row,row);
        
        for(i=0; i<=row-1; i++)
            tmp -= m_matLDLT(row,i) * m_matLDLT(i,row);
        
        if(tmp==(TT)0)
        {
            IM_THROW_NO_SOLUTION;
        }
        
        m_matLDLT(row,row) = tmp;
        
        for(i=row+1; i<size; i++)
        {
            TT tmp = m_matLDLT(i,row);
            
            for(k=0; k<row; k++)
                tmp -= m_matLDLT(i,k) * m_matLDLT(k,row);
            
            m_matLDLT(i,row) = tmp / m_matLDLT(row,row);
        }
    }
    
    core_block_copy_lower_tri(m_matLDLT.view().t(), m_matLDLT.view(), false);
}

// Solve Ax = y
template <typename TT>
im::Mtx<TT> im::MatrixDecompLDLT<TT>::solve(MtxView<TT> const &mavy) const
{
    IM_CHECK_VALID(mavy);
    IM_CHECK_ARGS(mavy.rows()==m_matLDLT.rows());
    
    int size = mavy.rows();
    int numcols = mavy.cols();
    
    Mtx<TT> mtxrtn(size, numcols);
    mtxrtn.view().copy_from(mavy);
    
    // backsub with matrix L (has 1s on diagonal)
    core_block_blas_trsm(mtxrtn.view(), m_matLDLT.view(), (TT)1, SideMode_L, TriMode_L, TransMode_N, DiagMode_U);
    
    // multiply by D^-1
    core_block_divide_rows(mtxrtn.view(), m_matLDLT.view().diag(), mtxrtn.view());
    
    // backsub with matrix L^T (has 1s on diagonal)
    core_block_blas_trsm(mtxrtn.view(), m_matLDLT.view(), (TT)1, SideMode_L, TriMode_U, TransMode_N, DiagMode_U);
    
    return mtxrtn;
}

template <typename TT>
im::Vec<TT> im::MatrixDecompLDLT<TT>::solve(VecView<TT> const &vvy) const
{
    IM_CHECK_VALID(vvy);
    IM_CHECK_ARGS(vvy.rows()==m_matLDLT.rows());
    
    int size = vvy.rows();
    
    Vec<TT> vrtn(size);
    vrtn.view().copy_from(vvy);
    
    // backsub with matrix L (has 1s on diagonal)
    core_block_blas_trsv(vrtn.view(), m_matLDLT.view(), TriMode_L, TransMode_N, DiagMode_U);
    
    // multiply by D^-1
    for(int i=0; i<size; i++)
        vrtn(i) /= m_matLDLT(i,i);
    
    // backsub with matrix L^T (has 1s on diagonal)
    core_block_blas_trsv(vrtn.view(), m_matLDLT.view(), TriMode_U, TransMode_N, DiagMode_U);
    
    return vrtn;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompLDLT<TT>::inverse() const
{
    int size = m_matLDLT.rows();
    
    Mtx<TT> matrtn(size, size);
    core_block_copy_lower_tri(matrtn.view(), m_matLDLT.view(), false);
    
    Vec<TT> vdiag(size);
    vdiag.copy_from(m_matLDLT.diag());
    
    matrtn.diag() = (TT)1;
    
    // Compute (L D L^T)^(-1) = L^(-T) D^(-1) L^(-1)
    
    // Invert lower triangle -> lower tri
    for(int i=size-2; i>=0; i--)
    {
        for(int j=size-1; j>=i+1; j--)
        {
            TT sum = (TT)0;
            for(int k=i+1; k<=j; k++)
                sum -= matrtn(j,k) * matrtn(k,i);
            matrtn(j,i) = sum;
        }
    }
    
    // Compute X^T D^(-1) X -> upper tri
    for(int i=0; i<size; i++)
    {
        for(int j=i; j<size; j++)
        {
            TT sum = (TT)0;
            for(int k=j; k<size; k++)
                sum +=  matrtn(k,i) * matrtn(k,j) / vdiag(k);
            matrtn(i,j) = sum;
        }
    }
    
    // Copy uppper to lower
    core_block_copy_lower_tri(matrtn.view(), matrtn.view().t(), false);
    
    return matrtn;
}

template <typename TT>
void im::MatrixDecompLLT<TT>::compute(MtxView<TT> const &mavA)
{
    int size = mavA.rows();
    
    IM_CHECK_VALID(mavA);
    IM_CHECK_MATRIX_SQUARE(mavA);
    IM_CHECK_ARGS(size>=1);
    
    m_matLLT.resize(size,size);
    
    core_block_copy_lower_tri(m_matLLT.view(), mavA, true);
    
    TT A00 = m_matLLT(0,0);
    if(A00<0)
    {
        IM_THROW_NO_SOLUTION;
    }
    TT L00 = std::sqrt(A00);
    m_matLLT(0,0) = L00;
    
    if(size>1)
    {
        TT L10 = m_matLLT(1,0) / L00;
        TT diag = m_matLLT(1,1) - L10 * L10;
        if(diag<0)
        {
            IM_THROW_NO_SOLUTION;
        }
        m_matLLT(1,0) = L10;
        m_matLLT(1,1) = std::sqrt(diag);
    }
    
    for(int k=2; k<size; k++)
    {
        for(int i=0; i<k; i++)
        {
            TT tmp = m_matLLT(k,i);
            if(i>0)
                tmp -= m_matLLT.block(i,0,1,i).dot_product(m_matLLT.block(k,0,1,i));
            
            m_matLLT(k,i) = tmp/m_matLLT(i,i);
        }
        
        TT diag = m_matLLT(k,k) - m_matLLT.block(k,0,1,k).magnitude_squared();
        if(diag<0)
        {
            IM_THROW_NO_SOLUTION;
        }
        m_matLLT(k,k) = std::sqrt(diag);
    }
    
    core_block_copy_lower_tri(m_matLLT.view().t(), m_matLLT.view(), false);
}

// Solve Ax = y
template <typename TT>
im::Mtx<TT> im::MatrixDecompLLT<TT>::solve(MtxView<TT> const &mavy) const
{
    int rows = mavy.rows();
    int cols = mavy.cols();
    
    IM_CHECK_VALID(mavy);
    IM_CHECK_ARGS(rows==m_matLLT.rows());
    
    Mtx<TT> matrtn(rows,cols);
    
    core_block_blas_trsm(matrtn.view(), m_matLLT.view(), (TT)1, SideMode_L, TriMode_L, TransMode_N, DiagMode_N);
    core_block_blas_trsm(matrtn.view(), m_matLLT.view(), (TT)1, SideMode_L, TriMode_U, TransMode_N, DiagMode_N);
    
    return matrtn;
}

template <typename TT>
im::Vec<TT> im::MatrixDecompLLT<TT>::solve(VecView<TT> const &vvy) const
{
    IM_CHECK_VALID(vvy);
    IM_CHECK_ARGS(vvy.rows()==m_matLLT.rows());
    
    int size = vvy.rows();
    
    Vec<TT> vrtn(size);
    vrtn.view().copy_from(vvy);
    
    core_block_blas_trsv(vrtn.view(), m_matLLT.view(), TriMode_L, TransMode_N, DiagMode_N);
    core_block_blas_trsv(vrtn.view(), m_matLLT.view(), TriMode_U, TransMode_N, DiagMode_N);
    
    return vrtn;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompLLT<TT>::inverse() const
{
    int size = m_matLLT.rows();
    
    Mtx<TT> matrtn(size, size);
    
    matrtn.copy_from(m_matLLT);
    
    // Compute (L L^T)^(-1) = L^(-T) L^(-1)
    
    // Invert lower triangle of LLT
    for(int j=size-1; j>=0; j--)
    {
        TT tmp = (TT)1 / matrtn(j,j);
        matrtn(j,j) = tmp;
        
        if(j<size-1)
        {
            int c = size - (j+1);
            VecView<TT> mv = matrtn.view().col(j).tail(c);
            core_block_blas_trmv(mv, matrtn.view().block(j+1, j+1, c, c), TriMode_L, TransMode_N, DiagMode_N);
            core_block_scale(mv, mv, -tmp);
        }
    }
    
    // Element ij is dot product of columns i and j of L^(-1) -> upper tri
    for(int i=0; i<size; i++)
    {
        for(int j = i+1; j<size; j++)
            matrtn(i,j) = matrtn.block(j,i,size-j,1).dot_product(matrtn.block(j,j,size-j,1));
        
        matrtn(i,i) = matrtn.block(i,i,size-i,1).magnitude_squared();
    }
    
    core_block_copy_lower_tri(matrtn.view(), matrtn.view().t(), false);
    
    return matrtn;
}

//TODO complex matrices

template class im::MatrixDecompLDLT<float>;
template class im::MatrixDecompLDLT<double>;

template class im::MatrixDecompLLT<float>;
template class im::MatrixDecompLLT<double>;

