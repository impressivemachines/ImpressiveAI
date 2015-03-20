//
//  block_loop.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/26/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_loop_h
#define Metaphor_block_loop_h

namespace im
{
    // Standard set of innner loops with 4 x unrolling.
    
    template <typename TT>
    inline void core_assign_loop(TT *pout, int strideout, TT const &x, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, pout += 4*strideout)
        {
            pout[0] = x;
            pout[strideout] = x;
            pout[2*strideout] = x;
            pout[3*strideout] = x;
        }
        for(; i<count; i++, pout += strideout)
            pout[0] = x;
    }
    
    template <typename TT>
    inline void core_copy_loop(TT *pout, int strideout, TT const *pin, int stridein, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, pout += 4*strideout, pin += 4*stridein)
        {
            pout[0] = pin[0];
            pout[strideout] = pin[stridein];
            pout[2*strideout] = pin[2*stridein];
            pout[3*strideout] = pin[3*stridein];
        }
        for(; i<count; i++, pout += strideout, pin += stridein)
            pout[0] = pin[0];
    }
    
    template <typename TT>
    inline void core_exchange_loop(TT *p1, int stride1, TT *p2, int stride2, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, p1 += 4*stride1, p2 += 4*stride2)
        {
            std::swap(p1[0], p2[0]);
            std::swap(p1[stride1], p2[stride2]);
            std::swap(p1[2*stride1], p2[2*stride2]);
            std::swap(p1[3*stride1], p2[3*stride2]);
        }
        for(; i<count; i++, p1 += stride1, p2 += stride2)
            std::swap(p1[0], p2[0]);
    }
    
    template <typename TT, typename TS>
    inline void core_madd_loop(TS &sum, TT const *p1, int stride1, TT const *p2, int stride2, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, p1 += 4*stride1, p2 += 4*stride2)
        {
            sum += p1[0] * p2[0];
            sum += p1[stride1] * p2[stride2];
            sum += p1[2*stride1] * p2[2*stride2];
            sum += p1[3*stride1] * p2[3*stride2];
        }
        for(; i<count; i++, p1 += stride1, p2 += stride2)
            sum += p1[0] * p2[0];
    }
    
    template <typename TT, typename TS>
    inline void core_maddconj_loop(TS &sum, TT const *p1, int stride1, TT const *p2, int stride2, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, p1 += 4*stride1, p2 += 4*stride2)
        {
            sum += core_conj(p1[0]) * p2[0];
            sum += core_conj(p1[stride1]) * p2[stride2];
            sum += core_conj(p1[2*stride1]) * p2[2*stride2];
            sum += core_conj(p1[3*stride1]) * p2[3*stride2];
        }
        for(; i<count; i++, p1 += stride1, p2 += stride2)
            sum += core_conj(p1[0]) * p2[0];
    }
    
    template <typename TT, typename TS>
    inline void core_msub_loop(TS &sum, TT const *p1, int stride1, TT const *p2, int stride2, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, p1 += 4*stride1, p2 += 4*stride2)
        {
            sum -= p1[0] * p2[0];
            sum -= p1[stride1] * p2[stride2];
            sum -= p1[2*stride1] * p2[2*stride2];
            sum -= p1[3*stride1] * p2[3*stride2];
        }
        for(; i<count; i++, p1 += stride1, p2 += stride2)
            sum -= p1[0] * p2[0];
    }
    
    template <typename TT, typename TS>
    inline void core_msubconj_loop(TS &sum, TT const *p1, int stride1, TT const *p2, int stride2, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, p1 += 4*stride1, p2 += 4*stride2)
        {
            sum -= core_conj(p1[0]) * p2[0];
            sum -= core_conj(p1[stride1]) * p2[stride2];
            sum -= core_conj(p1[2*stride1]) * p2[2*stride2];
            sum -= core_conj(p1[3*stride1]) * p2[3*stride2];
        }
        for(; i<count; i++, p1 += stride1, p2 += stride2)
            sum -= core_conj(p1[0]) * p2[0];
    }
    
    template <typename TT>
    inline void core_sadd_loop(TT *pout, int strideout, TT const *pin, int stridein, TT const &x, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, pout += 4*strideout, pin += 4*stridein)
        {
            pout[0] += x * pin[0];
            pout[strideout] += x * pin[stridein];
            pout[2*strideout] += x * pin[2*stridein];
            pout[3*strideout] += x * pin[3*stridein];
        }
        for(; i<count; i++, pout += strideout, pin += stridein)
            pout[0] += x * pin[0];
    }
    
    template <typename TT>
    inline void core_saddconj_loop(TT *pout, int strideout, TT const *pin, int stridein, TT const &x, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, pout += 4*strideout, pin += 4*stridein)
        {
            pout[0] += x * core_conj(pin[0]);
            pout[strideout] += x * core_conj(pin[stridein]);
            pout[2*strideout] += x * core_conj(pin[2*stridein]);
            pout[3*strideout] += x * core_conj(pin[3*stridein]);
        }
        for(; i<count; i++, pout += strideout, pin += stridein)
            pout[0] += x * core_conj(pin[0]);
    }
    
    template <typename TT>
    inline void core_ssub_loop(TT *pout, int strideout, TT const *pin, int stridein, TT const &x, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, pout += 4*strideout, pin += 4*stridein)
        {
            pout[0] -= x * pin[0];
            pout[strideout] -= x * pin[stridein];
            pout[2*strideout] -= x * pin[2*stridein];
            pout[3*strideout] -= x * pin[3*stridein];
        }
        for(; i<count; i++, pout += strideout, pin += stridein)
            pout[0] -= x * pin[0];
    }
    
    template <typename TT>
    inline void core_ssubconj_loop(TT *pout, int strideout, TT const *pin, int stridein, TT const &x, int count)
    {
        if(count<1)
            return;
        
        int const block4 = 4*(count/4);
        int i=0;
        for(; i<block4; i+=4, pout += 4*strideout, pin += 4*stridein)
        {
            pout[0] -= x * core_conj(pin[0]);
            pout[strideout] -= x * core_conj(pin[stridein]);
            pout[2*strideout] -= x * core_conj(pin[2*stridein]);
            pout[3*strideout] -= x * core_conj(pin[3*stridein]);
        }
        for(; i<count; i++, pout += strideout, pin += stridein)
            pout[0] -= x * core_conj(pin[0]);
    }
}

#endif
