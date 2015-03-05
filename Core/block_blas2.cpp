//
//  block_blas2.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/9/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT> void im::core_block_blas_gbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, int rowsA, int colsA, TT const &alpha, TT const &beta, int kl, int ku, TransMode trans)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_ARGS(mA.rows()>=rowsA);
    IM_CHECK_ARGS(mA.cols()>=1+ku+kl);
    
    int lenX, lenY, L, U;
    
    if(rowsA==0 || colsA==0)
        return;
    
    if(alpha == (TT)0 && beta == (TT)1)
        return;
    
    if(trans==TransMode_N)
    {
        lenX = colsA;
        lenY = rowsA;
        L = kl;
        U = ku;
    }
    else
    {
        lenX = rowsA;
        lenY = colsA;
        L = ku;
        U = kl;
    }
    
    IM_CHECK_ARGS(vx.rows()==lenX);
    IM_CHECK_ARGS(vy.rows()==lenY);
    
    if(beta==(TT)0)
    {
        vy = (TT)0;
        return;
    }
    else if(beta!=(TT)1)
    {
        for(int i=0; i<lenY; i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(trans==TransMode_N)
    {
        
        // y := alpha*A*x + y
        for(int i=0; i<lenY; i++)
        {
            TT temp = (TT)0;
            int j_min = std::max(i-L, 0);
            int j_max = std::min(i+U+1, lenX);
            
            if(j_max > j_min)
                core_madd_loop(temp, vx.ptr(j_min), vx.row_stride(), mA.ptr(i,L-i+j_min), mA.col_stride(), j_max - j_min);

            //for(int j=j_min; j<j_max; j++)
            //    temp += vx(j) * mA(i,L-i+j);
            
            vy(i) += alpha * temp;
        }
        
    }
    else if(trans==TransMode_T)
    {
        // y := alpha*A'*x + y
        for(int j=0; j<lenX; j++)
        {
            TT temp = alpha * vx(j);
            if(temp!=(TT)0)
            {
                int i_min = std::max(j-U, 0);
                int i_max = std::min(j+L+1, lenY);

                if(i_max > i_min)
                    core_sadd_loop(vy.ptr(i_min), vy.row_stride(), mA.ptr(j,U+i_min-j), mA.col_stride(), temp, i_max - i_min);
                
                //for(int i=i_min; i<i_max; i++)
                //    vy(i) += temp * mA(j,U+i-j);
            }
        }
    }
    else
    {
        // y := alpha*conj(A')*x + y
        for(int j=0; j<lenX; j++)
        {
            TT temp = alpha * vx(j);
            if(temp!=(TT)0)
            {
                int i_min = std::max(j-U, 0);
                int i_max = std::min(j+L+1, lenY);
                
                if(i_max > i_min)
                    core_saddconj_loop(vy.ptr(i_min), vy.row_stride(), mA.ptr(j,U+i_min-j), mA.col_stride(), temp, i_max - i_min);

                //for(int i=i_min; i<i_max; i++)
                //    vy(i) += temp * core_conj(mA(j,U+i-j));
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_gbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, int rowsA, int colsA, TT const &alpha, TT const &beta, int kl, int ku, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_gemv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TransMode trans)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    
    if(mA.rows()==0 || mA.cols()==0)
        return;
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    int lenX, lenY;
    
    if(trans==TransMode_N)
    {
        lenX = mA.cols();
        lenY = mA.rows();
    }
    else
    {
        lenX = mA.rows();
        lenY = mA.cols();
    }
    
    IM_CHECK_ARGS(vx.rows()==lenX);
    IM_CHECK_ARGS(vy.rows()==lenY);
    
    if(beta==(TT)0)
    {
        for(int i=0; i<lenY; i++)
            vy(i) = (TT)0;
    }
    else if(beta!=(TT)1)
    {
        for(int i=0; i<lenY; i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(trans==TransMode_N)
    {
        // y := alpha*A*x + y
        for(int i=0; i<lenY; i++)
        {
            TT temp = (TT)0;
            
            core_madd_loop(temp, mA.ptr(i,0), mA.col_stride(), vx.ptr(), vx.row_stride(), lenX);

            vy(i) += alpha * temp;
        }
    
    }
    else if(trans==TransMode_T)
    {
        // y := alpha*A'*x + y
        for(int i=0; i<lenY; i++)
        {
            TT temp = (TT)0;
            
            core_madd_loop(temp, mA.ptr(0,i), mA.row_stride(), vx.ptr(), vx.row_stride(), lenX);
            
            vy(i) += alpha * temp;
        }
    }
    else
    {
        // y := alpha*conj(A')*x + y
        for(int i=0; i<lenY; i++)
        {
            TT temp = (TT)0;
            
            core_maddconj_loop(temp, mA.ptr(0,i), mA.row_stride(), vx.ptr(), vx.row_stride(), lenX);
            
            //for(int j=0; j<lenX; j++)
            // temp += vx(j) * core_conj(mA(j,i));
            
            vy(i) += alpha * temp;
        }
    }
}

#define INST(TT) template void im::core_block_blas_gemv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_ger(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TransMode trans)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_ARGS(mA.rows()==vx.rows());
    IM_CHECK_ARGS(mA.cols()==vy.rows());
    
    if(alpha==(TT)0)
        return;
    
    if(trans==TransMode_C)
    {
        for(int j=0; j<mA.cols(); j++)
        {
            TT tmp = alpha * core_conj(vy(j));
            
            core_sadd_loop(mA.ptr(0,j), mA.row_stride(), vx.ptr(), vx.row_stride(), tmp, mA.rows());
            
            //for(int i=0; i<mA.rows(); i++)
            //    mA(i,j) += tmp * vx(i);
        }
    }
    else
    {
        for(int j=0; j<mA.cols(); j++)
        {
            TT tmp = alpha * vy(j);
            
            core_sadd_loop(mA.ptr(0,j), mA.row_stride(), vx.ptr(), vx.row_stride(), tmp, mA.rows());
            
            //for(int i=0; i<mA.rows(); i++)
            //    mA(i,j) += tmp * vx(i);
        }
    }
}

#define INST(TT) template void im::core_block_blas_ger(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_hbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    
    int sizeA = vx.rows();
    
    IM_CHECK_ARGS(mA.rows()>=sizeA);
    IM_CHECK_ARGS(mA.cols()>=1+k);
    IM_CHECK_ARGS(vy.rows()==sizeA);
    
    // y := alpha*A*x + beta*y
    if(sizeA==0)
        return;
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    if(beta==(TT)0)
        vy = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<vy.rows(); i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<sizeA; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            
            int j_min = i+1;
            int j_max = std::min(i+k+1, sizeA);
            
            TT Aii = mA(i,0);
            Aii.imag(0);
            vy(i) += temp1 * Aii;
            
            for(int j=j_min; j<j_max; j++)
            {
                TT Aij = mA(i,j-i);
                vy(j) += temp1 * core_conj(Aij);
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
    else
    {
        for(int i=0; i<sizeA; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            
            int j_min = std::max(i-k,0);
            int j_max = i;
            
            for(int j=j_min; j<j_max; j++)
            {
                TT Aij = mA(i,k-i+j);
                vy(j) += temp1 * core_conj(Aij);
                temp2 += vx(j) * Aij;
            }
            
            TT Aii = mA(i,0);
            Aii.imag(0);
            vy(i) += temp1 * Aii;
            vy(i) += alpha * temp2;
        }
    }
}

#define INST(TT) template void im::core_block_blas_hbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_hemv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_ARGS(mA.rows()==vy.rows());
    IM_CHECK_ARGS(mA.rows()==vx.rows());
    
    int N = mA.rows();
    
    // y := alpha*A*x + beta*y
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    if(beta==(TT)0)
        vy = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<N; i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            TT Aii = mA(i,0);
            Aii.imag(0);
            vy(i) += temp1 * Aii;
            
            for(int j=i+1; j<N; j++)
            {
                TT Aij = mA(i,j);
                vy(j) += temp1 * core_conj(Aij);
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            TT Aii = mA(i,0);
            Aii.imag(0);
            vy(i) += temp1 * Aii;
            
            for(int j=0; j<i; j++)
            {
                TT Aij = mA(i,j);
                vy(j) += temp1 * core_conj(Aij);
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
}

#define INST(TT) template void im::core_block_blas_hemv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_her(MtxView<TT> mA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_ARGS(mA.rows()==vx.rows());
    
    // A := alpha*x*conjg(x') + A
    
    if(alpha==(TT)0)
        return;
    
    int N = mA.rows();
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp = alpha * vx(i);
            
            mA(i,i) += temp * core_conj(vx(i));
            mA(i,i).imag(0); // ensure zero
            
            if(N-(i+1)>0)
                core_saddconj_loop(mA.ptr(i,i+1), mA.col_stride(), vx.ptr(i+1), vx.row_stride(), temp, N - (i+1));
            
            //for(int j=i+1; j<N; j++)
            //    mA(i,j) += temp * core_conj(vx(j));
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp = alpha * vx(i);
            
            core_saddconj_loop(mA.ptr(i,0), mA.col_stride(), vx.ptr(), vx.row_stride(), temp, i);
            
            //for(int j=0; j<i; j++)
            //    mA(i,j) += temp * core_conj(vx(j));
            
            mA(i,i) += temp * core_conj(vx(i));
            mA(i,i).imag(0); // ensure zero
        }
    }
    
}

#define INST(TT) template void im::core_block_blas_her(MtxView<TT> mA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_her2(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_ARGS(mA.rows()==vy.rows());
    IM_CHECK_ARGS(mA.rows()==vx.rows());
    
    // A := alpha * x*conjg(y') + conjg(alpha) * y*conjg(x') + A
    
    if(alpha==(TT)0)
        return;
    
    int N = mA.rows();
    
    TT conj_alpha = core_conj(alpha);

    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = conj_alpha * vy(i);

            mA(i,i).real(mA(i,i).real() + 2 * (temp1.real() * vy(i).real() + temp1.imag() * vy(i).imag()));
            mA(i,i).imag(0); // ensure zero
            
            for(int j=i+1; j<N; j++)
                mA(i,j) += temp1 * core_conj(vy(j)) + temp2 * core_conj(vx(j));
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = conj_alpha * vy(i);
            
            mA(i,i).real(mA(i,i).real() + 2 * (temp1.real() * vy(i).real() + temp1.imag() * vy(i).imag()));
            mA(i,i).imag(0); // ensure zero
            
            for(int j=0; j<i; j++)
                mA(i,j) += temp1 * core_conj(vy(j)) + temp2 * core_conj(vx(j));
        }
    }
}

#define INST(TT) template void im::core_block_blas_her2(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_hpmv(VecView<TT> vy, VecView<TT> const &vA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
   
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // y := alpha*A*x + beta*y
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    if(beta==(TT)0)
        vy = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<N; i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            TT Aii = vA(TPUP(N, i, i));
            Aii.imag(0);
            vy(i) += temp1 * Aii;
            
            for(int j=i+1; j<N; j++)
            {
                TT Aij = vA(TPUP(N, i, j));
                vy(j) += temp1 * core_conj(Aij);
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            TT Aii = vA(TPLO(N, i, i));
            Aii.imag(0);
            vy(i) += temp1 * Aii;
            
            for(int j=0; j<i; j++)
            {
                TT Aij = vA(TPLO(N, i, j));
                vy(j) += temp1 * core_conj(Aij);
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
}

#define INST(TT) template void im::core_block_blas_hpmv(VecView<TT> vy, VecView<TT> const &vA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_hpr(VecView<TT> vA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vx);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // A := alpha*x*conjg(x') + A
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp = alpha * vx(i);
            
            vA(TPUP(N, i, i)) += core_conj(vx(i)) * temp;
            vA(TPUP(N, i, i)).imag(0);
            
            for(int j=i+1; j<N; j++)
                vA(TPUP(N, i, j)) += core_conj(vx(j)) * temp;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp = alpha * vx(i);
            
            for(int j=0; j<i; j++)
                vA(TPLO(N, i, j)) += core_conj(vx(j)) * temp;
            
            vA(TPLO(N, i, i)) += core_conj(vx(i)) * temp;
            vA(TPLO(N, i, i)).imag(0);
        }
    }
}

#define INST(TT) template void im::core_block_blas_hpr(VecView<TT> vA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_hpr2(VecView<TT> vA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // A := alpha*x*conjg(y') + conjg(alpha)*y*conjg(x') + A
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = core_conj(alpha) * vy(i);
            
            vA(TPUP(N, i, i)).real(vA(TPUP(N, i, i)).real() + 2 * (temp1.real() * vy(i).real() + temp1.imag() * vy(i).imag()));
            vA(TPUP(N, i, i)).imag(0);
            
            for(int j=i+1; j<N; j++)
                vA(TPUP(N, i, j)) += temp1 * core_conj(vy(j)) + temp2 * core_conj(vx(j));
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = core_conj(alpha) * vy(i);
            
            vA(TPLO(N, i, i)).real(vA(TPLO(N, i, i)).real() + 2 * (temp1.real() * vy(i).real() + temp1.imag() * vy(i).imag()));
            vA(TPLO(N, i, i)).imag(0);
            
            for(int j=0; j<i; j++)
                vA(TPLO(N, i, j)) += temp1 * core_conj(vy(j)) + temp2 * core_conj(vx(j));
        }
    }
}

#define INST(TT) template void im::core_block_blas_hpr2(VecView<TT> vA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_sbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo)
{
    
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    
    int sizeA = vx.rows();
    
    IM_CHECK_ARGS(mA.rows()>=sizeA);
    IM_CHECK_ARGS(mA.cols()>=1+k);
    IM_CHECK_ARGS(vy.rows()==sizeA);
    
    // y := alpha*A*x + beta*y
    if(sizeA==0)
        return;
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    if(beta==(TT)0)
        vy = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<vy.rows(); i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<sizeA; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            
            int j_min = i+1;
            int j_max = std::min(i+k+1, sizeA);
            
            TT Aii = mA(i,0);
            vy(i) += temp1 * Aii;
            
            for(int j=j_min; j<j_max; j++)
            {
                TT Aij = mA(i,j-i);
                vy(j) += temp1 * Aij;
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
    else
    {
        for(int i=0; i<sizeA; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            
            int j_min = std::max(i-k,0);
            int j_max = i;
            
            for(int j=j_min; j<j_max; j++)
            {
                TT Aij = mA(i,k-i+j);
                vy(j) += temp1 * Aij;
                temp2 += vx(j) * Aij;
            }
            
            TT Aii = mA(i,0);
            vy(i) += temp1 * Aii;
            vy(i) += alpha * temp2;
        }
    }
}

#define INST(TT) template void im::core_block_blas_sbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_spmv(VecView<TT> vy, VecView<TT> const &vA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // y := alpha*A*x + beta*y
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    if(beta==(TT)0)
        vy = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<N; i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            TT Aii = vA(TPUP(N, i, i));
            vy(i) += temp1 * Aii;
            
            for(int j=i+1; j<N; j++)
            {
                TT Aij = vA(TPUP(N, i, j));
                vy(j) += temp1 * Aij;
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            TT Aii = vA(TPLO(N, i, i));
            vy(i) += temp1 * Aii;
            
            for(int j=0; j<i; j++)
            {
                TT Aij = vA(TPLO(N, i, j));
                vy(j) += temp1 * Aij;
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
}

#define INST(TT) template void im::core_block_blas_spmv(VecView<TT> vy, VecView<TT> const &vA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_spr(VecView<TT> vA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vx);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // A := alpha*x*x' + A
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp = alpha * vx(i);

            for(int j=i; j<N; j++)
                vA(TPUP(N, i, j)) += vx(j) * temp;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp = alpha * vx(i);
            
            for(int j=0; j<=i; j++)
                vA(TPLO(N, i, j)) += vx(j) * temp;
        }
    }
}

#define INST(TT) template void im::core_block_blas_spr(VecView<TT> vA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_spr2(VecView<TT> vA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // A := alpha*x*y' + conjg(alpha)*y*x' + A
    
    if(alpha==(TT)0)
        return;
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = alpha * vy(i);

            for(int j=i; j<N; j++)
                vA(TPUP(N, i, j)) += temp1 * vy(j) + temp2 * vx(j);
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = alpha * vy(i);

            for(int j=0; j<=i; j++)
                vA(TPLO(N, i, j)) += temp1 * vy(j) + temp2 * vx(j);
        }
    }
}

#define INST(TT) template void im::core_block_blas_spr2(VecView<TT> vA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_symv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    int N = mA.rows();
    
    IM_CHECK_ARGS(N==vx.rows());
    
    if(alpha==(TT)0 && beta==(TT)1)
        return;
    
    if(beta==(TT)0)
        vy = (TT)0;
    else if(beta!=(TT)1)
    {
        for(int i=0; i<N; i++)
            vy(i) *= beta;
    }
    
    if(alpha==(TT)0)
        return;
    
    // y := alpha*A*x + beta*y
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;
            
            vy(i) += temp1 * mA(i,i);
            
            for(int j=i+1; j<N; j++)
            {
                TT Aij = mA(i,j);
                vy(j) += temp1 * Aij;
                temp2 += vx(j) * Aij;
            }
            
            vy(i) += alpha * temp2;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT temp1 = alpha * vx(i);
            TT temp2 = (TT)0;

            vy(i) += temp1 * mA(i,i);
            
            for(int j=0; j<i; j++)
            {
                TT Aij = mA(i,j);
                vy(j) += temp1 * Aij;
                temp2 += vx(j) * Aij;
            }
        
            vy(i) += alpha * temp2;
        }
    }
}

#define INST(TT) template void im::core_block_blas_symv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_syr(MtxView<TT> mA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_ARGS(mA.rows()==vx.rows());
    
    if(alpha==(TT)0)
        return;
    
    // A := alpha*x*x' + A
    
    int N = mA.rows();
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT tmp = alpha * vx(i);
            
            core_sadd_loop(mA.ptr(i,i), mA.col_stride(), vx.ptr(i), vx.row_stride(), tmp, N-i);
            
            //for(int j=i; j<N; j++)
            //    mA(i,j) += vx(j) * tmp;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT tmp = alpha * vx(i);
            
            core_sadd_loop(mA.ptr(i,0), mA.col_stride(), vx.ptr(), vx.row_stride(), tmp, i+1);
            
            //for(int j=0; j<=i; j++)
            //    mA(i,j) += vx(j) * tmp;
        }
    }
}

#define INST(TT) template void im::core_block_blas_syr(MtxView<TT> mA, VecView<TT> const &vx, TT const &alpha, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_syr2(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vy);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_ARGS(mA.rows()==vx.rows());
    IM_CHECK_ARGS(mA.rows()==vy.rows());
    
    // A := alpha*x*y'+ alpha*y*x' + A
    
    int N = mA.rows();
    
    if(uplo==TriMode_U)
    {
        for(int i=0; i<N; i++)
        {
            TT tmp1 = alpha * vx(i);
            TT tmp2 = alpha * vy(i);
            
            for(int j=i; j<N; j++)
                mA(i,j) += vx(j) * tmp2 + vy(j) * tmp1;
        }
    }
    else
    {
        for(int i=0; i<N; i++)
        {
            TT tmp1 = alpha * vx(i);
            TT tmp2 = alpha * vy(i);
            
            for(int j=0; j<=i; j++)
                mA(i,j) += vx(j) * tmp2 + vy(j) * tmp1;
        }
    }
}

#define INST(TT) template void im::core_block_blas_syr2(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_tbmv(VecView<TT> vx, MtxView<TT> const &mA, DiagMode diag, int k, TriMode uplo, TransMode trans)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(mA.cols()>=N);
    IM_CHECK_ARGS(mA.rows()>=k+1);
    
    // x := A*x
    // or x:= A'*x
    // or x := conjg(A')*x
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT temp = (diag==DiagMode_N ? mA(i,0) : (TT)1.0) * vx(i);
                int j_min = i+1;
                int j_max = std::min(N,i+k+1);
                
                if(j_max > j_min)
                    core_madd_loop(temp, vx.ptr(j_min), vx.row_stride(), mA.ptr(i,j_min-i), mA.col_stride(), j_max - j_min);
                
                //for(int j=j_min; j<j_max; j++)
                //    temp += vx(j) * mA(i,j-i);
                
                vx(i) = temp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT temp = (diag==DiagMode_N ? mA(i,k) : (TT)1.0) * vx(i);
                int j_min = std::max(i-k,0);
                int j_max = i;
                
                if(j_max > j_min)
                    core_madd_loop(temp, vx.ptr(j_min), vx.row_stride(), mA.ptr(i,k-i+j_min), mA.col_stride(), j_max - j_min);
                
                //for(int j=j_min; j<j_max; j++)
                //    temp += vx(j) * mA(i,k-i+j);
                
                vx(i) = temp;
            }
        }
    }
    else if(trans==TransMode_T)
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT temp = (TT)0;
                int j_min = std::max(i-k,0);
                int j_max = i;
                
                for(int j=j_min; j<j_max; j++)
                    temp += vx(j) * mA(j,i-j);
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * mA(i,0);
                else
                    vx(i) += temp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT temp = (TT)0;
                int j_min = i+1;
                int j_max = std::min(N,i+k+1);
                
                for(int j=j_min; j<j_max; j++)
                    temp += vx(j) * mA(j,k-j+i);
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * mA(i,k);
                else
                    vx(i) += temp;
            }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT temp = (TT)0;
                int j_min = std::max(i-k,0);
                int j_max = i;
                
                for(int j=j_min; j<j_max; j++)
                    temp += vx(j) * core_conj(mA(j,i-j));
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * core_conj(mA(i,0));
                else
                    vx(i) += temp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT temp = (TT)0;
                int j_min = i+1;
                int j_max = std::min(N,i+k+1);
                
                for(int j=j_min; j<j_max; j++)
                    temp += vx(j) * core_conj(mA(j,k-j+i));
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * core_conj(mA(i,k));
                else
                    vx(i) += temp;
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_tbmv(VecView<TT> vx, MtxView<TT> const &mA, DiagMode diag, int k, TriMode uplo, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_tbsv(VecView<TT> vx, MtxView<TT> const &mA, DiagMode diag, int k, TriMode uplo, TransMode trans)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(mA.cols()>=N);
    IM_CHECK_ARGS(mA.rows()>=k+1);
    
    // A*x = b
    // or A'*x = b
    // or conjg(A')*x = b
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                int j_min = i+1;
                int j_max = std::min(N, i+k+1);
                
                if(j_max > j_min)
                    core_msub_loop(tmp, mA.ptr(i,j_min-i), mA.col_stride(), vx.ptr(j_min), vx.row_stride(), j_max - j_min);

                //for(int j=j_min; j<j_max; j++)
                //    tmp -= mA(i,j-i) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,0);
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
                int j_min = std::max(i-k,0);
                int j_max = i;
                
                if(j_max > j_min)
                    core_msub_loop(tmp, mA.ptr(i,k+j_min-i), mA.col_stride(), vx.ptr(j_min), vx.row_stride(), j_max - j_min);
                
                //for(int j=j_min; j<j_max; j++)
                //    tmp -= mA(i,k+j-i) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,k);
                else
                    vx(i) = tmp;
            }
        }
    }
    else if(trans==TransMode_T)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
                int j_min = std::max(i-k,0);
                int j_max = i;
                
                for(int j=j_min; j<j_max; j++)
                    tmp -= mA(j,i-j) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,0);
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                int j_min = i+1;
                int j_max = std::min(N, i+k+1);
                
                for(int j=j_min; j<j_max; j++)
                    tmp -= mA(j,k+i-j) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,k);
                else
                    vx(i) = tmp;
            }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
                int j_min = std::max(i-k,0);
                int j_max = i;
                
                for(int j=j_min; j<j_max; j++)
                    tmp -= core_conj(mA(j,i-j)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / core_conj(mA(i,0));
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                int j_min = i+1;
                int j_max = std::min(N, i+k+1);
                
                for(int j=j_min; j<j_max; j++)
                    tmp -= core_conj(mA(j,k+i-j)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / core_conj(mA(i,k));
                else
                    vx(i) = tmp;
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_tbsv(VecView<TT> vx, MtxView<TT> const &mA, DiagMode diag, int k, TriMode uplo, TransMode trans)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_tpmv(VecView<TT> vx, VecView<TT> const &vA, TriMode uplo, TransMode trans, DiagMode diag)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vx);
    
    int N = vx.rows();
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // x := A*x
    // or x := A'*x
    // or x := conjg(A')*x
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT atemp = vA(TPUP(N,i,i));
                TT temp = (diag==DiagMode_N ? atemp : (TT)1.0) * vx(i);
                
                for(int j=i+1; j<N; j++)
                    temp += vx(j) * vA(TPUP(N,i,j));
                
                vx(i) = temp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT atemp = vA(TPLO(N,i,i));
                TT temp = (diag==DiagMode_N ? atemp : (TT)1.0) * vx(i);

                for(int j=0; j<i; j++)
                    temp += vx(j) * vA(TPLO(N,i,j));
                
                vx(i) = temp;
            }
        }
    }
    else if(trans==TransMode_T)
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT atemp = vA(TPUP(N,i,i));
                TT temp = (diag==DiagMode_N ? atemp : (TT)1.0) * vx(i);
                
                for(int j=0; j<i; j++)
                    temp += vx(j) * vA(TPUP(N,j,i));
                
                vx(i) = temp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT atemp = vA(TPLO(N,i,i));
                TT temp = (diag==DiagMode_N ? atemp : (TT)1.0) * vx(i);
                
                for(int j=i+1; j<N; j++)
                    temp += vx(j) * vA(TPLO(N,j,i));

                vx(i) = temp;
            }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT atemp = core_conj(vA(TPUP(N,i,i)));
                TT temp = (diag==DiagMode_N ? atemp : (TT)1.0) * vx(i);
                
                for(int j=0; j<i; j++)
                    temp += vx(j) * core_conj(vA(TPUP(N,j,i)));
                
                vx(i) = temp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT atemp = core_conj(vA(TPLO(N,i,i)));
                TT temp = (diag==DiagMode_N ? atemp : (TT)1.0) * vx(i);
                
                for(int j=i+1; j<N; j++)
                    temp += vx(j) * core_conj(vA(TPLO(N,j,i)));
                
                vx(i) = temp;
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_tpmv(VecView<TT> vx, VecView<TT> const &vA, TriMode uplo, TransMode trans, DiagMode diag)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_tpsv(VecView<TT> vx, VecView<TT> const &vA, TriMode uplo, TransMode trans, DiagMode diag)
{
    IM_CHECK_VALID(vA);
    IM_CHECK_VALID(vx);
    
    int N = vx.rows();
    if(N<1)
        return;
    
    IM_CHECK_ARGS(vA.rows() >= N*(N-1)/2);
    
    // A*x = b
    // or A'*x = b
    // or conjg(A')*x = b
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);

                for(int j=i+1; j<N; j++)
                    tmp -= vA(TPUP(N,i,j)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / vA(TPUP(N,i,i));
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);

                for(int j=0; j<i; j++)
                    tmp -= vA(TPLO(N,i,j)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / vA(TPLO(N,i,i));
                else
                    vx(i) = tmp;
            }
        }
    }
    else if(trans==TransMode_T)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
                
                for(int j=0; j<i; j++)
                    tmp -= vA(TPUP(N,j,i)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / vA(TPUP(N,i,i));
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                
                for(int j=i+1; j<N; j++)
                    tmp -= vA(TPLO(N,j,i)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / vA(TPLO(N,i,i));
                else
                    vx(i) = tmp;
            }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
                
                for(int j=0; j<i; j++)
                    tmp -= core_conj(vA(TPUP(N,j,i))) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / core_conj(vA(TPUP(N,i,i)));
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                
                for(int j=i+1; j<N; j++)
                    tmp -= core_conj(vA(TPLO(N,j,i))) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / core_conj(vA(TPLO(N,i,i)));
                else
                    vx(i) = tmp;
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_tpsv(VecView<TT> vx, VecView<TT> const &vA, TriMode uplo, TransMode trans, DiagMode diag)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_trmv(VecView<TT> vx, MtxView<TT> const &mA, TriMode uplo, TransMode trans, DiagMode diag)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    int N = mA.rows();
    
    IM_CHECK_ARGS(vx.rows()==N);
    
    // x := A*x
    // or x := A'*x
    // or x := conjg(A')*x
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT temp = (TT)0;
                
                if(N-(i+1)>0)
                    core_madd_loop(temp, vx.ptr(i+1), vx.row_stride(), mA.ptr(i,i+1), mA.col_stride(), N - (i+1));
                
                //for(int j=i+1; j<N; j++)
                //    temp += vx(j) * mA(i,j);
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * mA(i,i);
                else
                    vx(i) += temp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT temp = (TT)0;
                
                core_madd_loop(temp, vx.ptr(), vx.row_stride(), mA.ptr(i,0), mA.col_stride(), i);
                
                //for(int j=0; j<i; j++)
                //    temp += vx(j) * mA(i,j);
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * mA(i,i);
                else
                    vx(i) += temp;
            }
        }
    }
    else if(trans==TransMode_T)
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT temp = (TT)0;
                
                core_madd_loop(temp, vx.ptr(), vx.row_stride(), mA.ptr(0,i), mA.row_stride(), i);
                
                //for(int j=0; j<i; j++)
                //    temp += vx(j) * mA(j,i);
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * mA(i,i);
                else
                    vx(i) += temp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT temp = (TT)0;
                
                if(N-(i+1)>0)
                    core_madd_loop(temp, vx.ptr(i+1), vx.row_stride(), mA.ptr(i+1,i), mA.row_stride(), N - (i+1));
                
                //for(int j=i+1; j<N; j++)
                //    temp += vx(j) * mA(j,i);
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * mA(i,i);
                else
                    vx(i) += temp;
            }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT temp = (TT)0;
                
                core_maddconj_loop(temp, mA.ptr(0,i), mA.row_stride(), vx.ptr(), vx.row_stride(), i);
                
                //for(int j=0; j<i; j++)
                //    temp += vx(j) * core_conj(mA(j,i));
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * core_conj(mA(i,i));
                else
                    vx(i) += temp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT temp = (TT)0;
                
                if(N-(i+1)>0)
                    core_maddconj_loop(temp, mA.ptr(i+1,i), mA.row_stride(), vx.ptr(i+1), vx.row_stride(), N-(i+1));
                
                //for(int j=i+1; j<N; j++)
                //    temp += vx(j) * core_conj(mA(j,i));
                
                if(diag==DiagMode_N)
                    vx(i) = temp + vx(i) * core_conj(mA(i,i));
                else
                    vx(i) += temp;
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_trmv(VecView<TT> vx, MtxView<TT> const &mA, TriMode uplo, TransMode trans, DiagMode diag)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_trsv(VecView<TT> vx, MtxView<TT> const &mA, TriMode uplo, TransMode trans, DiagMode diag)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    
    int N = mA.rows();
    
    IM_CHECK_ARGS(vx.rows()==N);
    
    // A*x = b
    // or A'*x = b
    // or conjg(A')*x = b
    
    if(trans==TransMode_N)
    {
        if(uplo==TriMode_U)
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                
                if(N-(i+1)>0)
                    core_msub_loop(tmp, mA.ptr(i,i+1), mA.col_stride(), vx.ptr(i+1), vx.row_stride(), N-(i+1));
                
                //for(int j=i+1; j<N; j++)
                //    tmp -= mA(i,j) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,i);
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
            
                core_msub_loop(tmp, mA.ptr(i,0), mA.col_stride(), vx.ptr(), vx.row_stride(), i);
                
                //for(int j=0; j<i; j++)
                //    tmp -= mA(i,j) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,i);
                else
                    vx(i) = tmp;
            }
        }
    }
    else if(trans==TransMode_T)
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);
                
                core_msub_loop(tmp, mA.ptr(0,i), mA.row_stride(), vx.ptr(), vx.row_stride(), i);
                
                //for(int j=0; j<i; j++)
                //    tmp -= mA(j,i) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,i);
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                
                if(N-(i+1)>0)
                    core_msub_loop(tmp, mA.ptr(i+1,i), mA.row_stride(), vx.ptr(i+1), vx.row_stride(), N - (i+1));

                //for(int j=i+1; j<N; j++)
                //    tmp -= mA(j,i) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / mA(i,i);
                else
                    vx(i) = tmp;
            }
        }
    }
    else
    {
        if(uplo==TriMode_U)
        {
            for(int i=0; i<N; i++)
            {
                TT tmp = vx(i);

                core_msubconj_loop(tmp, mA.ptr(0,i), mA.row_stride(), vx.ptr(), vx.row_stride(), i);

                //for(int j=0; j<i; j++)
                //    tmp -= core_conj(mA(j,i)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / core_conj(mA(i,i));
                else
                    vx(i) = tmp;
            }
        }
        else
        {
            for(int i=N-1; i>=0; i--)
            {
                TT tmp = vx(i);
                
                if(N-(i+1)>0)
                    core_msubconj_loop(tmp, mA.ptr(i+1,i), mA.row_stride(), vx.ptr(i+1), vx.row_stride(), N-(i+1));
                
                //for(int j=i+1; j<N; j++)
                //    tmp -= core_conj(mA(j,i)) * vx(j);
                
                if(diag==DiagMode_N)
                    vx(i) = tmp / core_conj(mA(i,i));
                else
                    vx(i) = tmp;
            }
        }
    }
}

#define INST(TT) template void im::core_block_blas_trsv(VecView<TT> vx, MtxView<TT> const &mA, TriMode uplo, TransMode trans, DiagMode diag)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

