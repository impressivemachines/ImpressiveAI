//
//  block_blas1.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/8/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_blas1_h
#define Metaphor_block_blas1_h

namespace im
{
    // Return sum of |xi|, or sum of |xi.re|+|xi.im| if complex
    // F, D, CF, CD
    template <typename TT> TT core_block_blas_asum(VecView<TT> const &vx);

    // y += alpha x
    // F, D, CF, CD
    template <typename TT> void core_block_blas_axpy(VecView<TT> vy, VecView<TT> const &vx, const TT &alpha);

    // y = x
    // F, D, CF, CD
    template <typename TT> void core_block_blas_copy(VecView<TT> vy, VecView<TT> const &vx);

    // Perform dot product (without conjugation, summation is equal to type)
    // F, D, CF, CD
    template <typename TT> TT core_block_blas_dot(VecView<TT> const &vx, VecView<TT> const &vy);
    
    // Perform dot product with conjugation of x
    // CF, CD
    template <typename TT> TT core_block_blas_dotc(VecView<TT> const &vx, VecView<TT> const &vy);
    
    // Double precision summation of dot product
    // F
    double core_block_blas_sdot(VecView<float> const &vx, VecView<float> const &vy);

    // Compute L2 norm
    // F, D, CF, CD
    template <typename TT> TT core_block_blas_nrm2(VecView<TT> const &vx);
    template <typename TT> std::complex<TT> core_block_blas_nrm2(VecView< std::complex<TT> > const &vx);
    template <typename TT> TT core_block_blas_nrm2(MtxView<TT> const &vx);
    template <typename TT> std::complex<TT> core_block_blas_nrm2(MtxView< std::complex<TT> > const &vx);
    
    // Performs rotation of points in the plane
    // xi = c*xi + s*yi
    // yi = c*yi - s*xi
    // F, D
    template <typename TT> void core_block_blas_rot(VecView<TT> vx, VecView<TT> vy, TT const &c, TT const &s);

    // Computes the parameters for a Givens rotation
    // F, D
    template <typename TT> void core_block_blas_rotg(TT &a, TT &b, TT &c, TT &s);

    // Performs modified Givens rotation of points in the plane
    // F, D
    template <typename TT> void core_block_blas_rotm(VecView<TT> vx, VecView<TT> vy, TT const *pparam);

    // Computes the parameters for a modified Givens rotation
    // F, D
    template <typename TT> void core_block_blas_rotmg(TT *pparam, TT &d1, TT &d2, TT &b1, TT const &b2);

    // Computes the product of a vector by a scalar
    // F, D, CF, CD
    template <typename TT> void core_block_blas_scal(VecView<TT> vx, TT const &alpha);

    // Swaps a vector with another vector
    // F, D, CF, CD
    template <typename TT> void core_block_blas_swap(VecView<TT> vx, VecView<TT> vy);
    
    // Finds the index of the element with maximum absolute value (complex abs is the sum of the abs of the re,im parts)
    // F, D, CF, CD
    template <typename TT> int core_block_blas_iamax(VecView<TT> const &vx);
    template <typename TT> int core_block_blas_iamax(VecView< std::complex<TT> > const &vx);
    
    // Finds the index of the element with minumum absolute value (complex abs is the sum of the abs of the re,im parts)
    // F, D, CF, CD
    template <typename TT> int core_block_blas_iamin(VecView<TT> const &vx);
    template <typename TT> int core_block_blas_iamin(VecView< std::complex<TT> > const &vx);
    
    // Computes absolute value of complex number |re| + |im| or a scalar |v|
    // F, D, CF, CD
    template <typename TT> TT core_block_blas_cabs1(std::complex<TT> const &v) { return core_sumabs(v); }
    
}
#endif
