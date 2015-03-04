//
//  block_blas3.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/9/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"


template <typename TT> void im::core_block_blas_gemm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, TransMode ta, TransMode tb)
{
    // Computes a matrix-matrix product with general matrices.
    // C := alpha*op(A)*op(B) + beta*C
    // op(X) is one of op(X) = X, or op(X) = XT, or op(X) = XH
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_VALID(mC);
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    MtxView<TT> mAA;
    if(ta==TransMode_N)
        mAA = mA;
    else
        mAA = mA.t();
    
    MtxView<TT> mBB;
    if(tb==TransMode_N)
        mBB = mB;
    else
        mBB = mB.t();
    
    IM_CHECK_ARGS(mC.rows()==mAA.rows());
    IM_CHECK_ARGS(mC.cols()==mBB.cols());
    IM_CHECK_ARGS(mAA.cols()==mBB.rows());
    
    if(beta==(TT)0)
    {
        mC = (TT)0;
    }
    else if(beta!=(TT)1)
    {
        for(int i=0; i<mC.rows(); i++)
            for(int j=0; j<mC.cols(); j++)
                mC(i,j) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(ta!=TransMode_C)
    {
        if(tb!=TransMode_C)
        {
            for(int k=0; k<mAA.cols(); k++)
                for(int i=0; i<mC.rows(); i++)
                {
                    TT temp = alpha * mAA(i,k);
                    if(temp!=(TT)0)
                    {
                        core_sadd_loop(mC.ptr(i,0), mC.col_stride(), mBB.ptr(k,0), mBB.col_stride(), temp, mC.cols());
                        
                        //for(int j=0; j<mC.cols(); j++)
                        //    mC(i,j) += temp * mBB(k,j);
                    }
                }
        }
        else
        {
            for(int k=0; k<mAA.cols(); k++)
                for(int i=0; i<mC.rows(); i++)
                {
                    TT temp = alpha * mAA(i,k);
                    if(temp!=(TT)0)
                    {
                        core_saddconj_loop(mC.ptr(i,0), mC.col_stride(), mBB.ptr(k,0), mBB.col_stride(), temp, mC.cols());
                        
                        //for(int j=0; j<mC.cols(); j++)
                        //    mC(i,j) += temp * core_conj(mBB(k,j));
                    }
                }
        }
    }
    else
    {
        if(tb!=TransMode_C)
        {
            for(int k=0; k<mAA.cols(); k++)
                for(int i=0; i<mC.rows(); i++)
                {
                    TT temp = alpha * core_conj(mAA(i,k));
                    if(temp!=(TT)0)
                    {
                        core_sadd_loop(mC.ptr(i,0), mC.col_stride(), mBB.ptr(k,0), mBB.col_stride(), temp, mC.cols());
                        
                        //for(int j=0; j<mC.cols(); j++)
                        //    mC(i,j) += temp * mBB(k,j);
                    }
                }
        }
        else
        {
            for(int k=0; k<mAA.cols(); k++)
                for(int i=0; i<mC.rows(); i++)
                {
                    TT temp = alpha * core_conj(mAA(i,k));
                    if(temp!=(TT)0)
                    {
                        core_saddconj_loop(mC.ptr(i,0), mC.col_stride(), mBB.ptr(k,0), mBB.col_stride(), temp, mC.cols());
                        
                        //for(int j=0; j<mC.cols(); j++)
                        //    mC(i,j) += temp * core_conj(mBB(k,j));
                    }
                }
        }
    }
}

#define INST(TT) template void im::core_block_blas_gemm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, TransMode ta, TransMode tb)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_blas_hemm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, SideMode side, TriMode uplo)
{
    // Computes a matrix-matrix product where input matrix A is Hermitian.
    // C := alpha*A*B + beta*C  (SideMode_L)
    // or
    // C := alpha*B*A + beta*C  (SideMode_R)
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_VALID(mC);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    int rows = mC.rows();
    int cols = mC.cols();
    
    if(beta==(TT)0)
        mC = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<rows; i++)
            for(int j=0; j<cols; j++)
                mC(i,j) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(side==SideMode_L)
    {
        // C += alpha*A*B
        IM_CHECK_ARGS(mA.rows()==rows);
        IM_CHECK_ARGS(mB.cols()==cols);
        IM_CHECK_ARGS(mA.cols()==mB.rows());

        if(uplo==TriMode_U)
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    
                    TT Aii = mA(i,i);
                    Aii.imag(0); // ensure zero
                    
                    mC(i,j) += temp1 * Aii;
                    
                    for(int k=i+1; k<rows; k++)
                    {
                        TT Aik = mA(i,k);
                        mC(k,j) += core_conj(Aik) * temp1;
                        temp2 += Aik * mB(k,j);
                    }
                    
                    mC(i,j) += alpha * temp2;
                }
        }
        else
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    
                    for(int k=0; k<i; k++)
                    {
                        TT Aik = mA(i,k);
                        mC(k,j) += core_conj(Aik) * temp1;
                        temp2 += Aik * mB(k,j);
                    }
                    
                    TT Aii = mA(i,i);
                    Aii.imag(0); // ensure zero
                    
                    mC(i,j) += temp1 * Aii;
                    mC(i,j) += alpha * temp2;
                }
        }
    }
    else
    {
        // C += alpha*B*A
        IM_CHECK_ARGS(mB.rows()==rows);
        IM_CHECK_ARGS(mA.cols()==cols);
        IM_CHECK_ARGS(mB.cols()==mA.rows());
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    
                    TT Ajj = mA(j,j);
                    Ajj.imag(0); // ensure zero
                    
                    mC(i,j) += temp1 * Ajj;
                    
                    for(int k=i+1; k<cols; k++)
                    {
                        TT Ajk = mA(j,k);
                        mC(i,k) += Ajk * temp1;
                        temp2 += core_conj(Ajk) * mB(i,k);
                    }
                    
                    mC(i,j) += alpha * temp2;
                }
        }
        else
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    
                    for(int k=0; k<j; k++)
                    {
                        TT Ajk = mA(j,k);
                        mC(i,k) += Ajk * temp1;
                        temp2 += core_conj(Ajk) * mB(i,k);
                    }
                    
                    TT Ajj = mA(j,j);
                    Ajj.imag(0); // ensure zero
                    
                    mC(i,j) += temp1 * Ajj;
                    mC(i,j) += alpha * temp2;
                }
        }
    }
}

#define INST(TT) template void im::core_block_blas_hemm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, SideMode side, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_herk(MtxView<TT> mC, MtxView<TT> const &mA, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
{
    // Performs a Hermitian rank-k update.
    // C := alpha*A*AH + beta*C (TransMode_N)
    // or
    // C := alpha*AH*A + beta*C (TransMode_C)
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mC);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_MATRIX_SQUARE(mC);
    IM_CHECK_ARGS(mA.rows()==mC.rows());
    
    if(beta==(TT)1 && (alpha==(TT)0 || kk==0))
        return;
    
    int N = mC.rows();
    
    if(beta==(TT)0)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) = (TT)0;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) = (TT)0;
        }
    }
    else if(beta!=(TT)1)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) *= beta;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) *= beta;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
            mC(i,i).imag(0);
    }
    
    if(alpha==(TT)0)
        return;
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp = (TT)0;
                    
                    core_maddconj_loop(temp, mA.ptr(j,0), mA.col_stride(), mA.ptr(i,0), mA.col_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += mA(i,k) * core_conj(mA(j,k));
                    
                    mC(i,j) += alpha * temp;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp = (TT)0;
                    
                    core_maddconj_loop(temp, mA.ptr(j,0), mA.col_stride(), mA.ptr(i,0), mA.col_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += mA(i,k) * core_conj(mA(j,k));
                    
                    mC(i,j) += alpha * temp;
                }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp = (TT)0;
                    
                    core_maddconj_loop(temp, mA.ptr(0,i), mA.row_stride(), mA.ptr(0,j), mA.row_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += core_conj(mA(k,i)) * mA(k,j);
                    
                    mC(i,j) += alpha * temp;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp = (TT)0;
                    
                    core_maddconj_loop(temp, mA.ptr(0,i), mA.row_stride(), mA.ptr(0,j), mA.row_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += core_conj(mA(k,i)) * mA(k,j);
                    
                    mC(i,j) += alpha * temp;
                }
        }
    }
}

#define INST(TT) template void im::core_block_blas_herk(MtxView<TT> mC, MtxView<TT> const &mA, TT const &alpha, TT const &beta, int k, TriMode uplo, TransMode trans)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_her2k(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
{
    // Performs a Hermitian rank-2k update.
    // C := alpha*A*BH + conjg(alpha)*B*AH + beta*C (TransMode_N)
    // or
    // C := alpha*BH*A + conjg(alpha)*AH*B + beta*C (TransMode_C)
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_VALID(mC);
    IM_CHECK_MATRIX_SQUARE(mC);
    
    if(beta==(TT)1 && (alpha==(TT)0 || kk==0))
        return;
    
    int N = mC.rows();
    
    if(beta==(TT)0)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) = (TT)0;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) = (TT)0;
        }
    }
    else if(beta!=(TT)1)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) *= beta;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) *= beta;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
            mC(i,i).imag(0);
    }
    
    if(alpha==(TT)0)
        return;
    
    TT conjalpha = core_conj(alpha);
    
    if(trans==TransMode_N)
    {
        IM_CHECK_ARGS(mA.rows()==N);
        IM_CHECK_ARGS(mA.cols()>=kk);
        IM_CHECK_ARGS(mB.rows()==N);
        IM_CHECK_ARGS(mB.cols()>=kk);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp1 = (TT)0;
                    TT temp2 = (TT)0;
                    for(int k=0; k<kk; k++)
                    {
                        temp1 += mA(i,k) * core_conj(mB(j,k));
                        temp2 += mB(i,k) * core_conj(mA(j,k));
                    }
                    mC(i,j) += alpha * temp1 + conjalpha * temp2;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp1 = (TT)0;
                    TT temp2 = (TT)0;
                    for(int k=0; k<kk; k++)
                    {
                        temp1 += mA(i,k) * core_conj(mB(j,k));
                        temp2 += mB(i,k) * core_conj(mA(j,k));
                    }
                    mC(i,j) += alpha * temp1 + conjalpha * temp2;
                }
        }
    }
    else
    {
        IM_CHECK_ARGS(mA.cols()==N);
        IM_CHECK_ARGS(mA.rows()>=kk);
        IM_CHECK_ARGS(mB.cols()==N);
        IM_CHECK_ARGS(mB.rows()>=kk);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp1 = (TT)0;
                    TT temp2 = (TT)0;
                    for(int k=0; k<kk; k++)
                    {
                        temp1 += core_conj(mA(k,i)) * mB(k,j);
                        temp2 += core_conj(mB(k,i)) * mA(k,j);
                    }
                    mC(i,j) += alpha * temp1 + conjalpha * temp2;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp1 = (TT)0;
                    TT temp2 = (TT)0;
                    for(int k=0; k<kk; k++)
                    {
                        temp1 += core_conj(mA(k,i)) * mB(k,j);
                        temp2 += core_conj(mB(k,i)) * mA(k,j);
                    }
                    mC(i,j) += alpha * temp1 + conjalpha * temp2;
                }
        }
    }
}

#define INST(TT) template void im::core_block_blas_her2k(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_symm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, SideMode side, TriMode uplo)
{
    // Computes a matrix-matrix product where input matrix A is symmetric.
    // C := alpha*A*B + beta*C  (SideMode_L)
    // or
    // C := alpha*B*A + beta*C  (SideMode_R)
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_VALID(mC);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    int rows = mC.rows();
    int cols = mC.cols();
    
    if(beta==(TT)0)
        mC = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<rows; i++)
            for(int j=0; j<cols; j++)
                mC(i,j) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(side==SideMode_L)
    {
        // AB
        IM_CHECK_ARGS(mA.rows()==rows);
        IM_CHECK_ARGS(mB.rows()==mA.cols());
        IM_CHECK_ARGS(mB.cols()==cols);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    mC(i,j) += temp1 * mA(i,i);
                    for(int k=i+1; k<rows; k++)
                    {
                        TT Aik = mA(i,k);
                        mC(k,j) += Aik * temp1;
                        temp2 += Aik * mB(k,j);
                    }
                    mC(i,j) += alpha * temp2;
                }
        }
        else
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    for(int k=0; k<i; k++)
                    {
                        TT Aik = mA(i,k);
                        mC(k,j) += Aik * temp1;
                        temp2 += Aik * mB(k,j);
                    }
                    mC(i,j) += temp1 * mA(i,i) + alpha * temp2;
                }
        }
    }
    else
    {
        // BA
        IM_CHECK_ARGS(mB.rows()==rows);
        IM_CHECK_ARGS(mB.cols()==mA.rows());
        IM_CHECK_ARGS(mA.cols()==cols);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    mC(i,j) += temp1 * mA(j,j);
                    for(int k=j+1; k<cols; k++)
                    {
                        TT Ajk = mA(j,k);
                        mC(i,k) += Ajk * temp1;
                        temp2 += Ajk * mB(i,k);
                    }
                    mC(i,j) += alpha * temp2;
                }
        }
        else
        {
            for(int i=0; i<rows; i++)
                for(int j=0; j<cols; j++)
                {
                    TT temp1 = alpha * mB(i,j);
                    TT temp2 = (TT)0;
                    for(int k=0; k<j; k++)
                    {
                        TT Ajk = mA(j,k);
                        mC(i,k) += Ajk * temp1;
                        temp2 += Ajk * mB(i,k);
                    }
                    mC(i,j) += temp1 * mA(j,j) + alpha * temp2;
                }
        }
    }
}

#define INST(TT) template void im::core_block_blas_symm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, SideMode side, TriMode uplo)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_syrk(MtxView<TT> mC, MtxView<TT> const &mA, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
{
    // Performs a symmetric rank-k update.
    // C := alpha*A*A' + beta*C (TransMode_N)
    // or
    // C := alpha*A'*A + beta*C (TransMode_T/C)
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mC);
    IM_CHECK_MATRIX_SQUARE(mC);
    
    if((alpha==(TT)0 || kk==0) && beta==(TT)1)
        return;
    
    int N = mC.rows();
    
    if(beta==(TT)0)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) = (TT)0;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) = (TT)0;
        }
    }
    else if(beta!=(TT)1)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) *= beta;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) *= beta;
        }
    }

    if(alpha==(TT)0)
        return;
    
    if(trans==TransMode_N)
    {
        IM_CHECK_ARGS(mA.cols()>=kk);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp = (TT)0;
                    
                    core_madd_loop(temp, mA.ptr(i,0), mA.col_stride(), mA.ptr(j,0), mA.col_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += mA(i,k) * mA(j,k);
                    mC(i,j) += alpha * temp;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp = (TT)0;
                    
                    core_madd_loop(temp, mA.ptr(i,0), mA.col_stride(), mA.ptr(j,0), mA.col_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += mA(i,k) * mA(j,k);
                    mC(i,j) += alpha * temp;
                }
        }
    }
    else
    {
        IM_CHECK_ARGS(mA.rows()>=kk);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp = (TT)0;
                    
                    core_madd_loop(temp, mA.ptr(0,i), mA.row_stride(), mA.ptr(0,j), mA.row_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += mA(k,i) * mA(k,j);
                    mC(i,j) += alpha * temp;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp = (TT)0;
                    
                    core_madd_loop(temp, mA.ptr(0,i), mA.row_stride(), mA.ptr(0,j), mA.row_stride(), kk);
                    
                    //for(int k=0; k<kk; k++)
                    //    temp += mA(k,i) * mA(k,j);
                    mC(i,j) += alpha * temp;
                }
        }
    }

}

#define INST(TT) template void im::core_block_blas_syrk(MtxView<TT> mC, MtxView<TT> const &mA, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_syr2k(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
{
    // Performs a symmetric rank-2k update.
    // C := alpha*A*B' + alpha*B*A' + beta*C (TransMode_N)
    // or
    // C := alpha*A'*B + alpha*B'*A + beta*C (TransMode_T/C)
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_VALID(mC);
    IM_CHECK_MATRIX_SQUARE(mC);
    
    if(beta==(TT)1 && (alpha==(TT)0 || kk==0))
        return;
    
    int N = mC.rows();
    
    if(beta==(TT)0)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) = (TT)0;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) = (TT)0;
        }
    }
    else if(beta!=(TT)1)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                    mC(i,j) *= beta;
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                    mC(i,j) *= beta;
        }
    }

    if(alpha==(TT)0)
        return;

    if(trans==TransMode_N)
    {
        IM_CHECK_ARGS(mA.rows()==N);
        IM_CHECK_ARGS(mA.cols()>=kk);
        IM_CHECK_ARGS(mB.rows()==N);
        IM_CHECK_ARGS(mB.cols()>=kk);
        
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
                for(int j=i; j<N; j++)
                {
                    TT temp = (TT)0;
                    for(int k=0; k<kk; k++)
                        temp += mA(i,k) * mB(j,k) + mB(i,k) * mA(j,k);
                    mC(i,j) += alpha * temp;
                }
        }
        else
        {
            for(int i=0; i<N; i++)
                for(int j=0; j<=i; j++)
                {
                    TT temp = (TT)0;
                    for(int k=0; k<kk; k++)
                        temp += mA(i,k) * mB(j,k) + mB(i,k) * mA(j,k);
                    mC(i,j) += alpha * temp;
                }
        }
    }
    else
    {
        IM_CHECK_ARGS(mA.cols()==N);
        IM_CHECK_ARGS(mA.rows()>=kk);
        IM_CHECK_ARGS(mB.cols()==N);
        IM_CHECK_ARGS(mB.rows()>=kk);
        
        if(uplo==TriMode_U)
        {
            for(int k=0; k<kk; k++)
                for(int i=0; i<N; i++)
                {
                    TT temp1 = alpha * mA(k,i);
                    TT temp2 = alpha * mB(k,i);
                    
                    for(int j=i; j<N; j++)
                        mC(i,j) += temp1 * mB(k,j) + temp2 * mA(k,j);
                }
        }
        else
        {
            for(int k=0; k<kk; k++)
                for(int i=0; i<N; i++)
                {
                    TT temp1 = alpha * mA(k,i);
                    TT temp2 = alpha * mB(k,i);
                    
                    for(int j=0; j<=i; j++)
                        mC(i,j) += temp1 * mB(k,j) + temp2 * mA(k,j);
                }
        }
    }
}

#define INST(TT) template void im::core_block_blas_syr2k(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_trmm(MtxView<TT> mB, MtxView<TT> const &mA, TT const &alpha, SideMode side, TriMode uplo, TransMode trans, DiagMode diag)
{
    // Computes a matrix-matrix product where one input matrix is triangular.
    // B := alpha*op(A)*B  (SideMode_L)
    // or
    // B := alpha*B*op(A)  (SideMode_R)
    // op(A) is one of op(A) = A, or op(A) = A', or op(A) = conjg(A').
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    MtxView<TT> mAA;
    if(trans==TransMode_N)
        mAA = mA;
    else
        mAA = mA.t();
    
    int rows = mB.rows();
    int cols = mB.cols();
    
    if(side==SideMode_L)
    {
        IM_CHECK_ARGS(mAA.rows()==rows);
        IM_CHECK_ARGS(mAA.cols()==rows);
        
        if(uplo==TriMode_U)
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * TriU(A) * B
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= mAA(i,i);
                        
                        if(rows-(i+1)>0)
                            core_madd_loop(temp, mAA.ptr(i,i+1), mAA.col_stride(), mB.ptr(i+1,j), mB.row_stride(), rows - (i+1));
                        
                        //for (int k = i + 1; k<rows; k++)
                        //    temp += mAA(i,k) * mB(k,j);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
            else
            {
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= core_conj(mAA(i,i));
                        
                        if(rows-(i+1)>0)
                            core_maddconj_loop(temp, mAA.ptr(i,i+1), mAA.col_stride(), mB.ptr(i+1,j), mB.row_stride(), rows - (i+1));
                        
                        //for (int k = i + 1; k<rows; k++)
                        //    temp += core_conj(mAA(i,k)) * mB(k,j);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
        }
        else
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * TriL(A) * B
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= mAA(i,i);
                        
                        core_madd_loop(temp, mAA.ptr(i,0), mAA.col_stride(), mB.ptr(0,j), mB.row_stride(), i);
                        
                       // for (int k = 0; k<i; k++)
                       //     temp += mAA(i,k) * mB(k,j);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
            else
            {
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= core_conj(mAA(i,i));
                        
                        core_maddconj_loop(temp, mAA.ptr(i,0), mAA.col_stride(), mB.ptr(0,j), mB.row_stride(), i);
                        
                        //for (int k = 0; k<i; k++)
                        //    temp += core_conj(mAA(i,k)) * mB(k,j);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
        }
    }
    else
    {
        IM_CHECK_ARGS(mAA.rows()==cols);
        IM_CHECK_ARGS(mAA.cols()==cols);
        
        if(uplo==TriMode_U)
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * B * TriU(A)
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= mAA(j,j);
                        
                        core_madd_loop(temp, mAA.ptr(0,j), mAA.row_stride(), mB.ptr(i,0), mB.col_stride(), j);
                        
                        //for (int k = 0; k<j; k++)
                        //    temp += mAA(k,j) * mB(i,k);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
            else
            {
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= core_conj(mAA(j,j));
                        
                        core_maddconj_loop(temp, mAA.ptr(0,j), mAA.row_stride(), mB.ptr(i,0), mB.col_stride(), j);
                        
                        //for (int k = 0; k<j; k++)
                        //    temp += core_conj(mAA(k,j)) * mB(i,k);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
        }
        else
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * B * TriL(A)
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= mAA(j,j);
                        
                        if(cols-(j+1)>0)
                            core_madd_loop(temp, mAA.ptr(j+1,j), mAA.row_stride(), mB.ptr(i,j+1), mB.col_stride(), cols - (j+1));
                        
                        //for (int k = j+1; k<cols; k++)
                        //    temp += mAA(k,j) * mB(i,k);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
            else
            {
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        TT temp = mB(i,j);
                        
                        if (diag==DiagMode_N)
                            temp *= core_conj(mAA(j,j));
                        
                        if(cols-(j+1)>0)
                            core_maddconj_loop(temp, mAA.ptr(j+1,j), mAA.row_stride(), mB.ptr(i,j+1), mB.col_stride(), cols - (j+1));
                        
                        //for (int k = j+1; k<cols; k++)
                        //    temp += core_conj(mAA(k,j)) * mB(i,k);
                        
                        mB(i,j) = alpha * temp;
                    }
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_trmm(MtxView<TT> mB, MtxView<TT> const &mA, TT const &alpha, SideMode side, TriMode uplo, TransMode trans, DiagMode diag)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_trsm(MtxView<TT> mB, MtxView<TT> const &mA, TT const &alpha, SideMode side, TriMode uplo, TransMode trans, DiagMode diag)
{
    // Solves a triangular matrix equation.
    // op(A)*X = alpha*B  (SideMode_L)
    // or
    // X*op(A) = alpha*B  (SideMode_R)
    // op(A) is one of op(A) = A, or op(A) = A', or op(A) = conjg(A').
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(mB);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    MtxView<TT> mAA;
    if(trans==TransMode_N)
        mAA = mA;
    else
        mAA = mA.t();
    
    int rows = mB.rows();
    int cols = mB.cols();
    
    if(alpha!=(TT)1)
    {
        for(int i=0; i<rows; i++)
            for(int j=0; j<cols; j++)
                mB(i,j) *= alpha;
    }
    
    if(side==SideMode_L)
    {
        IM_CHECK_ARGS(mAA.rows()==rows);
        IM_CHECK_ARGS(mAA.cols()==rows);
        
        if(uplo==TriMode_U)
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * inv(TriU(A)) * B
                for(int i=rows-1; i>=0; i--)
                {
                    if(diag==DiagMode_N)
                    {
                        TT Aii = mAA(i,i);
                        for(int j=0; j<cols; j++)
                            mB(i,j) /= Aii;
                    }
                    
                    for(int k=0; k<i; k++)
                    {
                        TT Aki = mAA(k,i);
                        
                        core_ssub_loop(mB.ptr(k,0), mB.col_stride(), mB.ptr(i,0), mB.col_stride(), Aki, cols);
                        
                        //for(int j=0; j<cols; j++)
                        //    mB(k,j) -= Aki * mB(i,j);
                    }
                }
            }
            else
            {
                for(int i=rows-1; i>=0; i--)
                {
                    if(diag==DiagMode_N)
                    {
                        TT Aii = core_conj(mAA(i,i));
                        for(int j=0; j<cols; j++)
                            mB(i,j) /= Aii;
                    }
                    
                    for(int k=0; k<i; k++)
                    {
                        TT Aki = core_conj(mAA(k,i));
                        
                        core_ssub_loop(mB.ptr(k,0), mB.col_stride(), mB.ptr(i,0), mB.col_stride(), Aki, cols);
                        
                        //for(int j=0; j<cols; j++)
                        //    mB(k,j) -= Aki * mB(i,j);
                    }
                }
            }
        }
        else
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * inv(TriL(A)) * B
                for(int i=0; i<rows; i++)
                {
                    if(diag==DiagMode_N)
                    {
                        TT Aii = mAA(i,i);
                        for(int j=0; j<cols; j++)
                            mB(i,j) /= Aii;
                    }
                    
                    for(int k=i+1; k<rows; k++)
                    {
                        TT Aki = mAA(k,i);
                        
                        core_ssub_loop(mB.ptr(k,0), mB.col_stride(), mB.ptr(i,0), mB.col_stride(), Aki, cols);
                        
                        //for(int j=0; j<cols; j++)
                        //    mB(k,j) -= Aki * mB(i,j);
                    }
                }
            }
            else
            {
                for(int i=0; i<rows; i++)
                {
                    if(diag==DiagMode_N)
                    {
                        TT Aii = core_conj(mAA(i,i));
                        for(int j=0; j<cols; j++)
                            mB(i,j) /= Aii;
                    }
                    
                    for(int k=i+1; k<rows; k++)
                    {
                        TT Aki = core_conj(mAA(k,i));
                        
                        core_ssub_loop(mB.ptr(k,0), mB.col_stride(), mB.ptr(i,0), mB.col_stride(), Aki, cols);
                        
                        //for(int j=0; j<cols; j++)
                        //    mB(k,j) -= Aki * mB(i,j);
                    }
                }
            }
        }
    }
    else
    {
        IM_CHECK_ARGS(mAA.rows()==cols);
        IM_CHECK_ARGS(mAA.cols()==cols);
        
        if(uplo==TriMode_U)
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * B * inv(TriU(A))
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        if(diag==DiagMode_N)
                            mB(i,j) /= mAA(j,j);
                        
                        TT Bij = mB(i,j);
                        
                        if(cols-(j+1)>0)
                            core_ssub_loop(mB.ptr(i,j+1), mB.col_stride(), mAA.ptr(j,j+1), mAA.col_stride(), Bij, cols - (j+1));
                        
                        //for(int k=j+1; k<cols; k++)
                        //    mB(i,k) -= mAA(j,k) * Bij;
                    }
            }
            else
            {
                for(int i=0; i<rows; i++)
                    for(int j=0; j<cols; j++)
                    {
                        if(diag==DiagMode_N)
                            mB(i,j) /= core_conj(mAA(j,j));
                        
                        TT Bij = mB(i,j);
                        
                        if(cols-(j+1)>0)
                            core_ssubconj_loop(mB.ptr(i,j+1), mB.col_stride(), mAA.ptr(j,j+1), mAA.col_stride(), Bij, cols - (j+1));
                        
                        //for(int k=j+1; k<cols; k++)
                        //    mB(i,k) -= core_conj(mAA(j,k)) * Bij;
                    }
            }
        }
        else
        {
            if(trans!=TransMode_C)
            {
                // B := alpha * B * inv(TriL(A))
                for(int i=0; i<rows; i++)
                    for(int j=cols-1; j>=0; j--)
                    {
                        if(diag==DiagMode_N)
                            mB(i,j) /= mAA(j,j);
                        
                        TT Bij = mB(i,j);
                        
                        core_ssub_loop(mB.ptr(i,0), mB.col_stride(), mAA.ptr(j,0), mAA.col_stride(), Bij, j);
                        
                        //for(int k=0; k<j; k++)
                        //    mB(i,k) -= mAA(j,k) * Bij;
                    }
            }
            else
            {
                for(int i=0; i<rows; i++)
                    for(int j=cols-1; j>=0; j--)
                    {
                        if(diag==DiagMode_N)
                            mB(i,j) /= core_conj(mAA(j,j));
                        
                        TT Bij = mB(i,j);
                        
                        core_ssubconj_loop(mB.ptr(i,0), mB.col_stride(), mAA.ptr(j,0), mAA.col_stride(), Bij, j);
                        
                        //for(int k=0; k<j; k++)
                        //    mB(i,k) -= core_conj(mAA(j,k)) * Bij;
                    }
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_trsm(MtxView<TT> mB, MtxView<TT> const &mA, TT const &alpha, SideMode side, TriMode uplo, TransMode trans, DiagMode diag)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

