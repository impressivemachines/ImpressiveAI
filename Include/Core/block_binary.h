//
//  block_binary.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_binary_h
#define Metaphor_block_binary_h

namespace im
{
    // All these operations are defined for:
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>

    // Binary operations (src can be the same as dst)

    
    // Element wise addition and subtraction multiply and divide
    template <typename TT> void core_block_add_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
    template <typename TT> void core_block_add_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    template <typename TT> void core_block_sub_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
    template <typename TT> void core_block_sub_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    template <typename TT> void core_block_multiply_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
    template <typename TT> void core_block_multiply_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    template <typename TT> void core_block_divide_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
    template <typename TT> void core_block_divide_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2);
    
    // Linearly blend between two vector element wise
    // dst = src1 if alpha = 0, src2 if alpha = 1
    template <typename TT> void core_block_blend(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2, float alpha);
    template <typename TT> void core_block_blend(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2, float alpha);
    
    // Add a scaled vector
    // dst = src1 + scale * src2
    template <typename TT> void core_block_add_scaled(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2, TT const &scale);
    template <typename TT> void core_block_add_scaled(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2, TT const &scale);

    // Misc
    
    // Outer product of two vectors
    template <typename TT> void core_block_outer_product(MtxView<TT> mdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
}

#endif


