//
//  block_blas3.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/9/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_blas3_h
#define Metaphor_block_blas3_h

namespace im
{
    // Matrix - Matrix products
    
    // Computes a matrix-matrix product with general matrices.
    // C := alpha*op(A)*op(B) + beta*C
    // op(X) is one of op(X) = X, or op(X) = XT, or op(X) = XH
    // F, D, CF, CD
    template <typename TT> void core_block_blas_gemm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, TransMode ta, TransMode tb);

    // Computes a matrix-matrix product where input matrix A is symmetric.
    // C := alpha*A*B + beta*C  (SideMode_L)
    // or
    // C := alpha*B*A + beta*C  (SideMode_R)
    // F, D, CF, CD
    template <typename TT> void core_block_blas_symm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, SideMode side, TriMode uplo);
    
    // Computes a matrix-matrix product where input matrix A is Hermitian.
    // C := alpha*A*B + beta*C  (SideMode_L)
    // or
    // C := alpha*B*A + beta*C  (SideMode_R)
    // CF, CD
    template <typename TT> void core_block_blas_hemm(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, SideMode side, TriMode uplo);
    
    // Computes a matrix-matrix product where one input matrix is triangular.
    // B := alpha*op(A)*B  (SideMode_L)
    // or
    // B := alpha*B*op(A)  (SideMode_R)
    // op(A) is one of op(A) = A, or op(A) = A', or op(A) = conjg(A').
    // F, D, CF, CD
    template <typename TT> void core_block_blas_trmm(MtxView<TT> mB, MtxView<TT> const &mA, TT const &alpha, SideMode side, TriMode uplo, TransMode trans, DiagMode diag);
    
    
    // Rank k updates
    
    // Performs a symmetric rank-k update.
    // C := alpha*A*A' + beta*C (TransMode_N)
    // or
    // C := alpha*A'*A + beta*C (TransMode_T/C)
    // F, D, CF, CD
    template <typename TT> void core_block_blas_syrk(MtxView<TT> mC, MtxView<TT> const &mA, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans);
    
    // Performs a Hermitian rank-k update.
    // C := alpha*A*AH + beta*C (TransMode_N)
    // or
    // C := alpha*AH*A + beta*C (TransMode_C)
    // CF, CD
    template <typename TT> void core_block_blas_herk(MtxView<TT> mC, MtxView<TT> const &mA, TT const &alpha, TT const &beta, int k, TriMode uplo, TransMode trans);
    
    
    // Rank 2k updates
    
    // Performs a symmetric rank-2k update.
    // C := alpha*A*B' + alpha*B*A' + beta*C (TransMode_N)
    // or
    // C := alpha*A'*B + alpha*B'*A + beta*C (TransMode_T/C)
    // F, D, CF, CD
    template <typename TT> void core_block_blas_syr2k(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans);
    
    // Performs a Hermitian rank-2k update.
    // C := alpha*A*BH + conjg(alpha)*B*AH + beta*C (TransMode_N)
    // or
    // C := alpha*BH*A + conjg(alpha)*AH*B + beta*C (TransMode_C)
    // CF, CD
    template <typename TT> void core_block_blas_her2k(MtxView<TT> mC, MtxView<TT> const &mA, MtxView<TT> const &mB, TT const &alpha, TT const &beta, int kk, TriMode uplo, TransMode trans);


    // Solve a triangular system
    
    // Solves a triangular matrix equation.
    // op(A)*X = alpha*B  (SideMode_L)
    // or
    // X*op(A) = alpha*B  (SideMode_R)
    // op(A) is one of op(A) = A, or op(A) = A', or op(A) = conjg(A').
    // F, D, CF, CD
    template <typename TT> void core_block_blas_trsm(MtxView<TT> mB, MtxView<TT> const &mA, TT const &alpha, SideMode side, TriMode uplo, TransMode trans, DiagMode diag);

}


#endif
