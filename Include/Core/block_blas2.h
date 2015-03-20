//
//  block_blas2.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/9/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_blas2_h
#define Metaphor_block_blas2_h

// Packed triangular matrix format
// N is number of rows/cols of square matrix
#define TRI_UP(N,i,j) (((i)*(2*(N)+1-(i)))/2 + (j) - (i))
#define TRI_LO(N,i,j) (((i)*((i)+1))/2 + (j))

namespace im
{
    // Matrix vector product using different types of matrix
    
    // Matrix-vector product using a general matrix
    // y := alpha*A*x + beta*y
    // or y := alpha*A*x + beta*y
    // or y := alpha*conjg(A')*x + beta*y
    // F, D, CF, CD
    template <typename TT> void core_block_blas_gemv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TransMode trans);

    // Matrix-vector product using a general band matrix
    // y := alpha*A*x + beta*y
    // or y := alpha*A'*x + beta*y
    // or y := alpha*conjg(A')*x + beta*y
    // F, D, CF, CD
    template <typename TT> void core_block_blas_gbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, int rowsA, int colsA, TT const &alpha, TT const &beta, int kl, int ku, TransMode trans);
    
    // Matrix-vector product using a symmetric matrix
    // y := alpha*A*x + beta*y
    // F, D
    template <typename TT> void core_block_blas_symv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo);
    
    // Matrix-vector product using symmetric band matrix
    // y := alpha*A*x + beta*y
    // F, D
    template <typename TT> void core_block_blas_sbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo);
    
    // Matrix-vector product using a symmetric *packed* matrix
    // y := alpha*A*x + beta*y,
    // F, D
    template <typename TT> void core_block_blas_spmv(VecView<TT> vy, VecView<TT> const &vA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo);

    // Matrix-vector product using a Hermitian matrix
    // y := alpha*A*x + beta*y
    // CF, CD
    template <typename TT> void core_block_blas_hemv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo);
    
    // Matrix-vector product using a Hermitian band matrix
    // y := alpha*A*x + beta*y
    // CF, CD
    template <typename TT> void core_block_blas_hbmv(VecView<TT> vy, MtxView<TT> const &mA, VecView<TT> const &vx, TT const &alpha, TT const &beta, int k, TriMode uplo);
    
    // Matrix-vector product using a Hermitian *packed* matrix
    // y := alpha*A*x + beta*y
    // CF, CD
    template <typename TT> void core_block_blas_hpmv(VecView<TT> vy, VecView<TT> const &vA, VecView<TT> const &vx, TT const &alpha, TT const &beta, TriMode uplo);

    // Matrix-vector product using a triangular matrix
    // x := A*x
    // or x := A'*x
    // or x := conjg(A')*x
    // F, D, CF, CD
    template <typename TT> void core_block_blas_trmv(VecView<TT> vx, MtxView<TT> const &mA, TriMode uplo, TransMode trans, DiagMode diag);
    
    // Matrix-vector product using a triangular band matrix
    // x := A*x
    // or x:= A'*x
    // or x := conjg(A')*x
    // F, D, CF, CD
    template <typename TT> void core_block_blas_tbmv(VecView<TT> vx, MtxView<TT> const &mA, DiagMode diag, int k, TriMode uplo, TransMode trans);
    
    // Matrix-vector product using a triangular *packed* matrix
    // x := A*x
    // or x := A'*x
    // or x := conjg(A')*x
    // F, D, CF, CD
    template <typename TT> void core_block_blas_tpmv(VecView<TT> vx, VecView<TT> const &vA, TriMode uplo, TransMode trans, DiagMode diag);
    
    
    
    // Rank 1 update of matrices

    // Rank-1 update of a general matrix
    // A := alpha*x*y' + A
    // or A := alpha*x*conjg(y') + A,
    // F, D, CF, CD
    template <typename TT> void core_block_blas_ger(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TransMode trans);

    // Rank-1 update of a symmetric matrix
    // A := alpha*x*x' + A
    // F, D
    template <typename TT> void core_block_blas_syr(MtxView<TT> mA, VecView<TT> const &vx, TT const &alpha, TriMode uplo);

    // Rank-1 update of a symmetric *packed* matrix
    // A:= alpha*x*x'+ A
    // F, D
    template <typename TT> void core_block_blas_spr(VecView<TT> vA, VecView<TT> const &vx, TT const &alpha, TriMode uplo);

    // Rank-1 update of a Hermitian matrix
    // A := alpha*x*conjg(x') + A
    // CF, CD
    template <typename TT> void core_block_blas_her(MtxView<TT> mA, VecView<TT> const &vx, TT const &alpha, TriMode uplo);

    // Rank-1 update of a Hermitian *packed* matrix
    // A := alpha*x*conjg(x') + A
    // CF, CD
    template <typename TT> void core_block_blas_hpr(VecView<TT> vA, VecView<TT> const &vx, TT const &alpha, TriMode uplo);
    
    
    
    // Rank2 update of matrices
    
    // Rank-2 update of a symmetric matrix
    // A := alpha*x*y'+ alpha*y*x' + A
    // F, D
    template <typename TT> void core_block_blas_syr2(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo);
    
    // Rank-2 update of a symmetric *packed* matrix
    // A:= alpha*x*y'+ alpha*y*x' + A
    // F, D
    template <typename TT> void core_block_blas_spr2(VecView<TT> vA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo);

    // Rank-2 update of a Hermitian matrix
    // A := alpha *x*conjg(y') + conjg(alpha)*y *conjg(x') + A
    // CF, CD
    template <typename TT> void core_block_blas_her2(MtxView<TT> mA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo);

    // Rank-2 update of a Hermitian *packed* matrix
    // A := alpha*x*conjg(y') + conjg(alpha)*y*conjg(x') + A
    // CF, CD
    template <typename TT> void core_block_blas_hpr2(VecView<TT> vA, VecView<TT> const &vx, VecView<TT> const &vy, TT const &alpha, TriMode uplo);



    // Solution of triangular linear system

    // Solution of a linear system of equations with a triangular matrix
    // A*x = b
    // or A'*x = b
    // or conjg(A')*x = b
    // F, D, CF, CD
    template <typename TT> void core_block_blas_trsv(VecView<TT> vx, MtxView<TT> const &mA, TriMode uplo, TransMode trans, DiagMode diag);

    // Solution of a linear system of equations with a triangular band matrix
    // A*x = b
    // or A'*x = b
    // or conjg(A')*x = b
    // F, D, CF, CD
    template <typename TT> void core_block_blas_tbsv(VecView<TT> vx, MtxView<TT> const &mA, DiagMode diag, int k, TriMode uplo, TransMode trans);

    // Solution of a linear system of equations with a triangular *packed* matrix
    // A*x = b
    // or A'*x = b
    // or conjg(A')*x = b
    // F, D, CF, CD
    template <typename TT> void core_block_blas_tpsv(VecView<TT> vx, VecView<TT> const &vA, TriMode uplo, TransMode trans, DiagMode diag);





}

#endif
