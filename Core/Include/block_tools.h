//
//  block_tools.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/25/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_tools_h
#define Metaphor_block_tools_h

namespace im
{

    // Copy/clear actions on upper/lower triangular part of matrix
    // The source/destination matrices do not need to be square or of the same size
    // They only have to be big enough to hold the clear/copy region (including the diagonal)
    template <typename TT> void core_block_clear_upper_tri(MtxView<TT> mdst, bool include_diagonal);
    template <typename TT> void core_block_clear_lower_tri(MtxView<TT> mdst, bool include_diagonal);
    template <typename TT> void core_block_copy_upper_tri(MtxView<TT> mdst, MtxView<TT> const &msrc, bool include_diagonal);
    template <typename TT> void core_block_copy_lower_tri(MtxView<TT> mdst, MtxView<TT> const &msrc, bool include_diagonal);
    
    
    
    // Normalize by computing the sum of the squares of elements of src and computing dst = src / (k+std::sqrt(sum_squares))
    template <typename TT> void core_block_normalize(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &k = (TT)0);
    template <typename TT> void core_block_normalize(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &k = (TT)0);
    
    // Independently normalize all columns (use transpose to normlize rows)
    template <typename TT> void core_block_normalize_columns(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &k = (TT)0);
    
    // Compute the matrix vector product xT A x
    template <typename TT> TT core_block_compute_xTAx(MtxView<TT> const &mA, VecView<TT> const &vx);
    
    // vScale is a vector of scale factors
    // Multiplies every row by corresponding elements in vScale(row)
    template <typename TT> void core_block_multiply_rows(MtxView<TT> mdst, VecView<TT> const &vScale, MtxView<TT> const &msrc);
    
    // vScale is a vector of scale factors
    // Divides every row by corresponding elements in vScale(row)
    template <typename TT> void core_block_divide_rows(MtxView<TT> mdst, VecView<TT> const &vScale, MtxView<TT> const &msrc);
    
    // Selection of best match column
    // Returns the column index for the column vector in mSet which has the closest L2 distance to mVector
    template <typename TT> int core_block_nearest(MtxView<TT> const &mSet, VecView<TT> const &v);
}

#endif
