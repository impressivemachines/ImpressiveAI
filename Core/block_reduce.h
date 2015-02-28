//
//  block_reduce.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_block_reduce_h
#define ImpressiveAI_block_reduce_h

namespace im
{
    // Reduction operations that distill a matrix down to a single value

    // All these operations are defined for:
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>
    
    // Vectors
    
    // Add all elements together
    template <typename TT> TT core_block_reduce_add(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_add(MtxView<TT> const &msrc);
    
    // Multiply all elements
    template <typename TT> TT core_block_reduce_multiply(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_multiply(MtxView<TT> const &msrc);
    
    // Compute the mean of all the elements
    template <typename TT> TT core_block_reduce_mean(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_mean(MtxView<TT> const &msrc);
    
    // Compute the variance and optionally the mean of all elements
    template <typename TT> TT core_block_reduce_variance(VecView<TT> const &vsrc, TT *pmeanrtn = NULL);
    template <typename TT> TT core_block_reduce_variance(MtxView<TT> const &msrc, TT *pmeanrtn = NULL);
    
    // Compute the sum of the squares of all the elements
    template <typename TT> TT core_block_reduce_sum_squares(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_sum_squares(MtxView<TT> const &msrc);
    
    // Compute the sum of the squares of the differences between the vector elements
    template <typename TT> TT core_block_reduce_squared_distance(VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
    template <typename TT> TT core_block_reduce_squared_distance(MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    // Compute the dot product between the vectors by multiplying corresponding elements and forming the sum
    template <typename TT> TT core_block_reduce_multiply_add(VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
    template <typename TT> TT core_block_reduce_multiply_add(MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);


    // The reduce operations below work with
    //   uint8_t
    //   int16_t
    //   int
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>

    template <typename TT> TT core_block_reduce_max(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_max(MtxView<TT> const &msrc);
    
    template <typename TT> TT core_block_reduce_min(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_min(MtxView<TT> const &msrc);
    
    template <typename TT> TT core_block_reduce_max_abs(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_max_abs(MtxView<TT> const &msrc);
    
    template <typename TT> TT core_block_reduce_min_abs(VecView<TT> const &vsrc);
    template <typename TT> TT core_block_reduce_min_abs(MtxView<TT> const &msrc);
    
    // Locate the element with rank order k. If k=-1 then the median is returned
    // use k = 0 for smallest, k = size-1 for largest
    template <typename TT> TT core_block_reduce_median(VecView<TT> const &vsrc, int k = -1);
    template <typename TT> TT core_block_reduce_median(MtxView<TT> const &msrc, int k = -1);
    
    
    // Reductions within rows, by combining data across columns, to create a vector with one column and as many rows as the source matrix
    // To reduce over down columns, merely transpose the source matrix before calling the function

    // The row reduce operations below work with
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>
    
    // Add all the elements of each row independently
    template <typename TT> void core_block_reduce_rows_add(VecView<TT> vdst, MtxView<TT> const &msrc);
    
    // Form the product of all elements on each row independenly
    template <typename TT> void core_block_reduce_rows_multiply(VecView<TT> vdst, MtxView<TT> const &msrc);
    
    // Compute the mean of each row independently
    template <typename TT> void core_block_reduce_rows_mean(VecView<TT> vdst, MtxView<TT> const &msrc);
    
    // Compute the variance of each row independently
    template <typename TT> void core_block_reduce_rows_variance(VecView<TT> vdst, MtxView<TT> const &msrc);
    
    // Compute the sum of the squares of all the each row independently
    template <typename TT> void core_block_reduce_rows_sum_squares(VecView<TT> vdst, MtxView<TT> const &msrc);
    
    // Compute the sum of the squares of the differences between the matrix elements for each row independently
    template <typename TT> void core_block_reduce_rows_squared_distance(VecView<TT> vdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    // Dot produc of each row between the two matrices
    template <typename TT> void core_block_reduce_rows_multiply_add(VecView<TT> vdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    
    // The row reduce operations below work with
    //   uint8_t
    //   int16_t
    //   int
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>
    
    template <typename TT> void core_block_reduce_rows_max(VecView<TT> vdst, MtxView<TT> const &msrc);
    template <typename TT> void core_block_reduce_rows_min(VecView<TT> vdst, MtxView<TT> const &msrc);
    template <typename TT> void core_block_reduce_rows_max_abs(VecView<TT> vdst, MtxView<TT> const &msrc);
    template <typename TT> void core_block_reduce_rows_min_abs(VecView<TT> vdst, MtxView<TT> const &msrc);
    
    // Locate the element with rank order k. If k=-1 then the median is returned
    template <typename TT> void core_block_reduce_rows_median(VecView<TT> vdst, MtxView<TT> const &msrc, int k = -1);
    
}

#endif

