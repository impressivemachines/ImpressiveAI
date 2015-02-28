//
//  decomp_qr.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::MatrixDecompQR<TT>::compute(MtxView<TT> const &mavA, bool use_pivoting)
{
    m_unpacked = false;
    m_unpacked_thin = false;
    
    IM_CHECK_VALID(mavA);
    
    int M = mavA.rows();
    int N = mavA.cols();
    
    int i,j;
    
    m_matQ.resize(0,0); // Will be M x M
    m_matR.resize(0,0); // Will be M x N
    
    m_matQR.resize(M, N);
    m_matQR.copy_from(mavA);
    
    Vec<TT> vNorm;
    if(use_pivoting)
        vNorm.resize(N);
    
    m_perm.resize(N);
    
    int size = std::min(M,N);
    m_vTau.resize(size);
    
    if(use_pivoting)
        for(i=0; i<N; i++)
            vNorm(i) = std::sqrt(core_block_reduce_sum_squares(m_matQR.view().col(i)));
    
    for(i=0; i<size; i++)
    {
        if(use_pivoting)
        {
            int loc = vNorm.block(i,N-i).location_of_max();
            loc += i;
            
            if(loc != i)
            {
                core_block_exchange(m_matQR.view().col(i), m_matQR.view().col(loc));
                m_perm.swap(i, loc);
                std::swap(vNorm(i), vNorm(loc));
            }
        }
        
        VecView<TT> vHH = m_matQR.view().col(i).block(i,M-i); // column down from i,i
        
        TT tau = core_compute_householder_vector(vHH);
        m_vTau(i) = tau;
        
        if(i+1 < N)
            core_apply_householder_vector(m_matQR.view().block(i,i+1,M-i,N-(i+1)), tau, vHH);
        
        if(use_pivoting)
        {
            if(i+1<M)
            {
                for(j=i+1; j<N; j++)
                {
                    TT x = vNorm(j);
                    if(x > (TT)0)
                    {
                        TT temp = m_matQR(i,j) / x;
                        
                        TT y;
                        if(std::abs(temp)>=(TT)1)
                            y = (TT)0;
                        else
                            y = x * std::sqrt((TT)1-temp*temp);
                        
                        // recompute the norm
                        if(std::abs(y/x) < std::numeric_limits<TT>::epsilon() * (TT)4.472) // 4.472 = sqrt(20)
                            y = m_matQR.block(i+1,j,M-(i+1),1).magnitude();
                        
                        vNorm(j) = y;
                    }
                }
            }
        }
    }
}

template <typename TT>
void im::MatrixDecompQR<TT>::unpack_QR()
{
    if(m_unpacked)
        return;
    
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    int size = std::min(M,N);
    
    m_matQ.resize(M, M);
    m_matR.resize(M, N);
    
    m_matQ.set_identity();
    
    int i,j;
    
    // Compute Q
    for(i=size-1; i>=0; i--)
        core_apply_householder_vector(m_matQ.view().block(i,i,M-i,M-i), m_vTau(i), m_matQR.view().col(i).block(i,M-i));
    
    // Compute R
    for(i=0; i<M; i++)
    {
        for(j=0; j<i && j<N; j++)
            m_matR(i,j) = (TT)0;
        for(j=i; j<N; j++)
            m_matR(i,j) = m_matQR(i,j);
    }
    
    m_unpacked = true;
    m_unpacked_thin = false;
}

template <typename TT>
void im::MatrixDecompQR<TT>::unpack_QR_thin()
{
    if(m_unpacked_thin)
        return;
    
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    int size = std::min(M,N);
    
    m_matQ.resize(M, size);
    m_matR.resize(size, N);
    
    m_matQ.set_identity();
    
    int i,j;
    
    // Compute Q
    for(i=size-1; i>=0; i--)
        core_apply_householder_vector(m_matQ.view().block(i,i,M-i,size-i), m_vTau(i), m_matQR.view().col(i).block(i,M-i));
    
    // Compute R
    for(i=0; i<size; i++)
    {
        for(j=0; j<i && j<N; j++)
            m_matR(i,j) = (TT)0;
        for(j=i; j<N; j++)
            m_matR(i,j) = m_matQR(i,j);
    }
    
    m_unpacked_thin = true;
    m_unpacked = false;
}

// Solve Ax = y for matrix x and y using QRPT factorization
template <typename TT>
im::Mtx<TT> im::MatrixDecompQR<TT>::solve(MtxView<TT> const &mavy) const
{
    IM_CHECK_VALID(mavy);
    
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    
    IM_CHECK_ARGS(M==N);
    IM_CHECK_ARGS(M==mavy.rows());
    
    Mtx<TT> matx(M, mavy.cols());
    matx.copy_from(mavy);
    
    for(int col=0; col<mavy.cols(); col++)
    {
        // Compute s = Q^T y
        for(int i=0; i<M; i++)
            core_apply_householder_vector(matx.view().block(i,col,M-i,1), m_vTau(i), m_matQR.view().col(i).block(i,M-i));
        
        // Solve R b = s
        core_block_blas_trsv(matx.view().col(col), m_matQR.view(), TriMode_U, TransMode_N, DiagMode_N);
        
        // Solve P^T x = b by inverse permutation
        m_perm.matrix_PTX_in_place(matx.view().col(col).matrix_view());
    }
    
    return matx;
}

template <typename TT>
im::Vec<TT> im::MatrixDecompQR<TT>::solve(VecView<TT> const &vvy) const
{
    IM_CHECK_VALID(vvy);
    
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    
    IM_CHECK_ARGS(M==N);
    IM_CHECK_ARGS(M==vvy.rows());
    
    Vec<TT> vx(M);
    vx.copy_from(vvy);
    
    // Compute s = Q^T y
    for(int i=0; i<M; i++)
        core_apply_householder_vector(vx.view().block(i,M-i).matrix_view(), m_vTau(i), m_matQR.view().col(i).block(i,M-i));
    
    // Solve R b = s
    core_block_blas_trsv(vx.view(), m_matQR.view(), TriMode_U, TransMode_N, DiagMode_N);
    
    // Solve P^T x = b by inverse permutation
    m_perm.matrix_PTX_in_place(vx.view().matrix_view());
    
    return vx;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompQR<TT>::solve_least_squares(MtxView<TT> const &mavy) const
{
    IM_CHECK_VALID(mavy);
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    int cols = mavy.cols();
    IM_CHECK_ARGS(mavy.rows()==M);
    IM_CHECK_ARGS(M>=N);
    int size = std::min(M,N);
    
    Mtx<TT> matx(N,cols);
    Vec<TT> vtmp(M);
    
    for(int col=0; col<cols; col++)
    {
        vtmp.copy_from(mavy.col(col));
        
        // Compute s = Q^T y
        for(int i=0; i<size; i++)
            core_apply_householder_vector(vtmp.view().block(i,M-i).matrix_view(), m_vTau(i), m_matQR.view().col(i).block(i,M-i));
        
        VecView<TT> vtop = vtmp.view().head(N);
        
        // Solve R b = s
        core_block_blas_trsv(vtop, m_matQR.view().block(0,0,N,N), TriMode_U, TransMode_N, DiagMode_N);
        
        // Solve P^T x = b by inverse permutation
        m_perm.matrix_PTX_in_place(vtop.matrix_view());
        
        matx.view().col(col).copy_from(vtop);
    }
    
    return matx;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompQR<TT>::inverse() const
{
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    
    Mtx<TT> matI(M, N);
    if(M!=N)
    {
        matI = (TT)0;
        return matI;
    }
    
    matI.set_identity();
    
    return solve(matI.view());
}

template <typename TT>
void im::MatrixDecompQR<TT>::matrix_QTX_in_place(MtxView<TT> mavX) const
{
    IM_CHECK_VALID(mavX);
    int M = m_matQR.rows();
    int N = m_matQR.cols();
    int size = std::min(M,N);
    int cols = mavX.cols();
    IM_CHECK_ARGS(mavX.rows()==M);
    
    for(int i=0; i<size; i++)
        core_apply_householder_vector(mavX.block(i,0,M-i,cols), m_vTau(i), m_matQR.view().col(i).block(i,M-i));
}

template class im::MatrixDecompQR<float>;
template class im::MatrixDecompQR<double>;

