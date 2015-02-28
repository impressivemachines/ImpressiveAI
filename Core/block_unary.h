//
//  block_unary.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_block_unary_h
#define ImpressiveAI_block_unary_h

namespace im
{
    // All these operations are defined for:
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>

    // Vector/Matrix View unary operations
    
    template <typename TT> void core_block_neg(VecView<TT> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_block_neg(MtxView<TT> mdst, MtxView<TT> const &msrc);
    
    template <typename TT> void core_block_conj(VecView<TT> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_block_conj(MtxView<TT> mdst, MtxView<TT> const &msrc);
    
    template <typename TT> void core_block_recip(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &epsilon = TT(0)); // 1/(x+epsilon)
    template <typename TT> void core_block_recip(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &epsilon = TT(0)); // 1/(x+epsilon)
    
    template <typename TT> void core_block_abs(VecView<TT> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_block_abs(MtxView<TT> mdst, MtxView<TT> const &msrc);
    
    template <typename TT> void core_block_sqrt(VecView<TT> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_block_sqrt(MtxView<TT> mdst, MtxView<TT> const &msrc);
    
    template <typename TT> void core_block_pow(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &exponent);
    template <typename TT> void core_block_pow(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &exponent);
    
    template <typename TT> void core_block_exp(VecView<TT> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_block_exp(MtxView<TT> mdst, MtxView<TT> const &msrc);
    
    template <typename TT> void core_block_log(VecView<TT> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_block_log(MtxView<TT> mdst, MtxView<TT> const &msrc);
    
    template <typename TT> void core_block_sigm(VecView<TT> vdst, VecView<TT> const &vsrc); // y = 1/(1+std::exp(-x))
    template <typename TT> void core_block_sigm(MtxView<TT> mdst, MtxView<TT> const &msrc); // y = 1/(1+std::exp(-x))
    
    template <typename TT> void core_block_scale(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale);
    template <typename TT> void core_block_scale(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale);
    
    template <typename TT> void core_block_offset(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &offet);
    template <typename TT> void core_block_offset(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &offet);
    
    template <typename TT> void core_block_scale_offset(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale, TT const &offset); // y = x * scale + offset
    template <typename TT> void core_block_scale_offset(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale, TT const &offset); // y = x * scale + offset

}

#endif
