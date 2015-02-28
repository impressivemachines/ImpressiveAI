//
//  utils.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_utils_h
#define ImpressiveAI_utils_h

namespace im
{
        
    // Handle endian-ness for file io
    // Only works for integer types
    template <typename T> T core_byte_reverse(T w)
    {
        T r = 0;
        for(int i=0; i<sizeof(r); i++)
        {
            r<<=8;
            r |= w & 0xff;
            w >>=8;
        }
        return r;
    }

    template <typename T> void core_byte_reverse_in_place(T &w)
    {
        w = core_byte_reverse(w);
    }

    inline bool core_is_big_endian()
    {
        union { uint64_t quad; uint32_t islittle; } t;
        t.quad = 1;
        
        return t.islittle==0;
    }

    template <typename T> T core_to_big_endian(T w)
    {
        if(core_is_big_endian())
            return w;
        else
            return core_byte_reverse(w);
    }

    template <typename T> T core_to_little_endian(T w)
    {
        if(core_is_big_endian())
            return core_byte_reverse(w);
        else
            return w;
    }

    // Special case sorting for n variables
    // Available for int, float, double
    template <typename TT> void core_sort_2(TT &a, TT &b, SortDirection dir);
    template <typename TT> void core_sort_3(TT &a, TT &b, TT &c, SortDirection dir);
    template <typename TT> void core_sort_4(TT &a, TT &b, TT &c, TT &d, SortDirection dir);

    // Test for equality within defined precision
    template <typename TT> bool core_equal(TT x, TT y, double epsilon)
    {
        return std::abs(x-y)<=epsilon;
    }

    template <typename TT> bool core_equal(VecView<TT> const &vv1, VecView<TT> const &vv2, double epsilon)
    {
        IM_CHECK_VALID(vv1);
        IM_CHECK_VALID(vv2);
        IM_CHECK_VECTOR_SIZES_MATCH(vv1, vv2);
        
        for(int r=0; r<vv1.rows(); r++)
            if(!im::core_equal(vv1(r), vv2(r), epsilon))
                return false;
        return true;
    }
    
    template <typename TT> bool core_equal(MtxView<TT> const &mav1, MtxView<TT> const &mav2, double epsilon)
    {
        IM_CHECK_VALID(mav1);
        IM_CHECK_VALID(mav2);
        IM_CHECK_MATRIX_SIZES_MATCH(mav1, mav2);
        
        for(int r=0; r<mav1.rows(); r++)
            for(int c=0; c<mav1.cols(); c++)
                if(!im::core_equal(mav1(r,c), mav2(r,c), epsilon))
                    return false;
        return true;
    }
    
    // Ensure that an angle lies in the range -PI < angle <= PI
    template <typename T> T core_cyclic_modulus(T angle);

    // Modulus that handles negative numbers for val, wrapping into the positive range
    inline int core_modulus(int val, int mod)
    {
        if(val>=0)
            return val % mod;
        else
            return val + mod * ((mod - val - 1)/mod);
    }
    
    // clamp to range 0 - (mod-1)
    inline int core_clamp(int val, int mod)
    {
        val = std::max(val, 0);
        val = std::min(val, mod-1);
        return val;
    }
    
    // triangular reflection of range, e.g. with mod=4:
    // in = -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    // out = 2,  1, 0, 1, 2, 3, 2, 1, 0, 1, 2, 3
    // out always lies from 0 through mod-1
    inline int core_reflect(int val, int mod)
    {
        if(mod<2)
            return 0;
        val = core_modulus(val, 2*mod-2);
        if(val>=mod)
            val = 2*(mod-1) - val;
        return val;
    }
    
    
    // Calculate log to the base of 2. Returns -1 if value is not a power of 2.
    int core_integer_log2(uint64_t value);

    // Returns true if vale is a power of 2
    inline bool core_is_pow2(uint64_t value) { return (value & (value-1))==0; }

    // Compute m = sqrt(x*x+y*y) with better precision
    template <typename TT> TT core_hypot(TT x, TT y);

    // get absolute time in seconds - source in time.cpp
    double core_get_time();
}

#endif
