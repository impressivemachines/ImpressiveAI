//
//  decomp_svd.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::MatrixDecompSVD<TT>::in_place_jacobi(MtxView<TT> mvA, MtxView<TT> mvQ, VecView<TT> vvS) const
{
    int M = mvA.rows();
    int N = mvA.cols();
    
    IM_CHECK_ARGS(M>=N);
    IM_CHECK_MATRIX_SQUARE(mvQ);
    IM_CHECK_ARGS(mvQ.rows()==N);
    IM_CHECK_VECTOR_SIZE(vvS, N);
    
    mvQ = (TT)0;
    mvQ.diag() = (TT)1;
    
    // store column error estimates in S for use during orthogonalization
    for(int j=0; j<N; j++)
        vvS(j) = TypeProperties<TT>::epsilon() * core_block_blas_nrm2(mvA.col(j));
    
    int count = 1;
    int sweep = 0;
    int sweepmax = std::max(5*N, 12);
    
    TT tolerance = 10 * M * TypeProperties<TT>::epsilon();
    
    // orthogonalize A by plane rotations
    while(count>0 && sweep <= sweepmax)
    {
        //printf("%d\n", sweep);
        
        count = N * (N-1) / 2;
        
        for(int j=0; j<N-1; j++)
        {
            for(int k=j+1; k<N; k++)
            {
                VecView<TT> cj = mvA.col(j);
                VecView<TT> ck = mvA.col(k);
                
                TT p = (TT)2 * core_block_blas_dot(cj, ck);
                TT a = core_block_blas_nrm2(cj);
                TT b = core_block_blas_nrm2(ck);
                TT q = a*a-b*b;
                TT v = core_hypot(p, q);
                
                // test for cols j,k orthogonal, or dominant errors
                
                bool sorted = (core_coerce(a) >= core_coerce(b));
                bool orthog = (std::abs(p) <= tolerance * core_coerce(a*b));
                
                TT abserr_a = vvS(j);
                TT abserr_b = vvS(k);
                
                if(sorted && (orthog || a < abserr_a || b < abserr_b))
                {
                    count--;
                    continue;
                }
                
                TT cosine, sine;
                
                // compute rotation angles
                if(v==(TT)0 || !sorted)
                {
                    cosine = (TT)0;
                    sine = (TT)1;
                }
                else
                {
                    cosine = std::sqrt((v+q) / ((TT)2 * v));
                    sine = p / ((TT)2 * v * cosine);
                }
                
                // apply rotation to A
                for(int i=0; i<M; i++)
                {
                    TT Aik = mvA(i,k);
                    TT Aij = mvA(i,j);
                    mvA(i,j) = Aij * cosine + Aik * sine;
                    mvA(i,k) = -Aij * sine + Aik * cosine;
                }
                
                vvS(j) = std::abs(cosine) * abserr_a + std::abs(sine) * abserr_b;
                vvS(k) = std::abs(sine) * abserr_a + std::abs(cosine) * abserr_b;
                
                // apply rotation to Q
                for(int i=0; i<N; i++)
                {
                    TT Qij = mvQ(i,j);
                    TT Qik = mvQ(i,k);
                    mvQ(i,j) = Qij * cosine + Qik * sine;
                    mvQ(i,k) = -Qij * sine + Qik * cosine;
                }
            }
        }
        
        sweep++;
    }
    
    // compute singular values
    TT prev_norm = (TT)(-1);
    
    for(int j=0; j<N; j++)
    {
        VecView<TT> vcol = mvA.col(j);
        TT norm = core_block_blas_nrm2(vcol);
        
        // test if singular value is zero
        if(norm==(TT)0 || prev_norm==(TT)0 || (j>0 && norm <= tolerance * prev_norm))
        {
            vvS(j) = (TT)0;
            vcol = (TT)0;
            prev_norm = (TT)0;
        }
        else
        {
            vvS(j) = norm;
            core_block_scale(vcol, vcol, (TT)1 / norm);
            prev_norm = norm;
        }
    }
    
    if(count>0)
        IM_THROW_NO_SOLUTION;
}

template <typename TT>
void im::MatrixDecompSVD<TT>::compute(MtxView<TT> const &mvA, bool qr_precondition)
{
    IM_CHECK_VALID(mvA);
    
    int M = mvA.rows();
    int N = mvA.cols();
    
    IM_CHECK_ARGS(M>=N);
    
    m_mU.resize(M,N);
    m_mV.resize(N,N);
    m_vS.resize(N);
    
    m_mU.view().copy_from(mvA);
    
    if(N==1)
    {
        // svd of a column vector
        VecView<TT> vcol = m_mU.view().col(0);
        TT norm = core_block_blas_nrm2(vcol);
        m_vS(0) = norm;
        m_mV(0,0) = (TT)1;
        if(norm!=(TT)0)
            core_block_scale(vcol, vcol, (TT)1 / norm);
        return;
    }
    
    if(qr_precondition)
    {
        // precondition using QR decomposition:
        // A = QR
        // SVD of R: A = Q (USV')
        // A = (QU) S V'
        
        // convert A into upper triangular matrix R
        for(int i=0; i<N; i++)
        {
            VecView<TT> vHH = m_mU.view().col(i).block(i,M-i); // column down from i,i
            
            TT tau = core_compute_householder_vector(vHH);
            m_vS(i) = tau;
            
            if(i+1 < N)
                core_apply_householder_vector(m_mU.view().block(i,i+1,M-i,N-(i+1)), tau, vHH);
        }
        
        Mtx<TT> mX(N,N);
        
        // copy upper tri part of U into X
        for(int i=0; i<N; i++)
        {
            int j;
            for(j=0; j<i; j++)
                mX(i,j) = (TT)0;
            for(; j<N; j++)
                mX(i,j) = m_mU(i,j);
        }
        
        // Convert U into orthogonal matrix L
        for(int j=N-1; j >= 0; j--)
            core_apply_householder_identity(m_mU.view().block(j,j,M-j,N-j), m_vS(j));
        
        // SVD of X
        in_place_jacobi(mX.view(), m_mV.view(), m_vS.view());
        
        // obtain U = LX
        Vec<TT> vsum(N);
        
        for(int i=0; i<M; i++)
        {
            vsum = (TT)0;
            
            for(int j=0; j<N; j++)
                core_block_blas_axpy(vsum.view(), mX.view().row(j), m_mU(i,j));
            
            m_mU.view().row(i).copy_from(vsum.view());
        }
    }
    else
    {
        in_place_jacobi(m_mU.view(), m_mV.view(), m_vS.view());
    }
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompSVD<TT>::solve(MtxView<TT> const &mavy) const
{
    IM_CHECK_VALID(mavy);
    IM_CHECK_ARGS(mavy.rows()==m_mU.rows());
    
    int N = m_mU.cols();
    
    Mtx<TT> mrtn(N, mavy.cols());
    Vec<TT> vtmp(N);
    
    for(int c=0; c<mavy.cols(); c++)
    {
        core_block_blas_gemv(vtmp.view(), m_mU.view(), mavy.col(c), (TT)1, (TT)0, TransMode_T);
        
        for(int i=0; i<N; i++)
        {
            TT alpha = m_vS(i);
            if(alpha!=(TT)0)
                alpha = (TT)1/alpha;
            vtmp(i) *= alpha;
        }
        
        core_block_blas_gemv(mrtn.view().col(c), m_mV.view(), vtmp.view(), (TT)1, (TT)0, TransMode_N);
    }
    
    return mrtn;
}

template <typename TT>
im::Vec<TT> im::MatrixDecompSVD<TT>::solve(VecView<TT> const &vvy) const
{
    IM_CHECK_VALID(vvy);
    IM_CHECK_ARGS(vvy.rows()==m_mU.rows());
    
    int N = m_mU.cols();
    
    Vec<TT> vtmp(N);
    Vec<TT> vrtn(N);
    
    core_block_blas_gemv(vtmp.view(), m_mU.view(), vvy, (TT)1, (TT)0, TransMode_T);
    
    for(int i=0; i<N; i++)
    {
        TT alpha = m_vS(i);
        if(alpha!=(TT)0)
            alpha = (TT)1/alpha;
        vtmp(i) *= alpha;
    }
    
    core_block_blas_gemv(vrtn.view(), m_mV.view(), vtmp.view(), (TT)1, (TT)0, TransMode_N);
    
    return vrtn;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompSVD<TT>::pseudo_inverse() const
{
    int M = m_mU.rows();
    int N = m_mU.cols();
    
    Mtx<TT> mrtn(N, M);
    Vec<TT> vtmp(N);
    
    for(int c=0; c<M; c++)
    {
        vtmp.view().copy_from(m_mU.view().row(c));

        for(int i=0; i<N; i++)
        {
            TT alpha = m_vS(i);
            if(alpha!=(TT)0)
                alpha = (TT)1/alpha;
            vtmp(i) *= alpha;
        }
        
        core_block_blas_gemv(mrtn.view().col(c), m_mV.view(), vtmp.view(), (TT)1, (TT)0, TransMode_N);
    }

    return mrtn;
}

template class im::MatrixDecompSVD<float>;
template class im::MatrixDecompSVD<double>;


template <typename TT>
void im::core_decomp_svd_2x2(MtxView<TT> mvU, VecView<TT> vvS, MtxView<TT> mvV, MtxView<TT> const &mvA)
{
    IM_CHECK_VALID(mvU);
    IM_CHECK_VALID(vvS);
    IM_CHECK_VALID(mvV);
    IM_CHECK_VALID(mvA);
    IM_CHECK_VECTOR_SIZE(vvS, 2);
    IM_CHECK_MATRIX_SIZE(mvA, 2, 2);
    IM_CHECK_MATRIX_SIZE(mvU, 2, 2);
    IM_CHECK_MATRIX_SIZE(mvV, 2, 2);
    
    TT data[4];
    MtxView<TT> mAAT(2,2,2,1,data);

    core_block_blas_gemm(mAAT, mvA, mvA, (TT)1, (TT)0, TransMode_N, TransMode_T);
    
    core_decomp_eigen_2x2(vvS, mvU, mAAT);
    
    vvS(0) = std::sqrt(vvS(0));
    vvS(1) = std::sqrt(vvS(1));
    
    // V = (S^-1 U^T A)^T
    core_block_blas_gemm(mvV.t(), mvU, mvA, (TT)1, (TT)0, TransMode_T, TransMode_N);
    
    if(vvS(0)>(TT)0)
    {
        mvV(0,0) /= vvS(0);
        mvV(1,0) /= vvS(0);
    }
    else
    {
        mvV(0,0) = (TT)1;
        mvV(1,0) = (TT)0;
    }
    
    if(vvS(1)>(TT)0)
    {
        mvV(0,1) /= vvS(1);
        mvV(1,1) /= vvS(1);
    }
    else
    {
        mvV(0,1) = -mvV(1,0);
        mvV(1,1) = mvV(0,0);
    }
}

template void im::core_decomp_svd_2x2(MtxView<float> mvU, VecView<float> vvS, MtxView<float> mvV, MtxView<float> const &mvA);
template void im::core_decomp_svd_2x2(MtxView<double> mvU, VecView<double> vvS, MtxView<double> mvV, MtxView<double> const &mvA);


