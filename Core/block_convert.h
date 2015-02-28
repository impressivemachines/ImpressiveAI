//
//  block_convert.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/25/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_block_convert_h
#define ImpressiveAI_block_convert_h

namespace im
{
    // Cast conversion inline functions for templates
    inline void core_type_convert(float &dst, double const &src) { dst = (float)src; }
    inline void core_type_convert(double &dst, float const &src) { dst = (double)src; }
    inline void core_type_convert(Cf &dst, Cd const &src) { dst = (Cf)src; }
    inline void core_type_convert(Cd &dst, Cf const &src) { dst = (Cd)src; }
    inline void core_type_convert(Cf &dst, float const &src) { dst = (Cf)src; }
    inline void core_type_convert(Cd &dst, double const &src) { dst = (Cd)src; }
    inline void core_type_convert(float &dst, Cf const &src) { dst = src.real(); } // just takes real part
    inline void core_type_convert(double &dst, Cd const &src) { dst = src.real(); } // just takes real part

    // Conversion from any integer type
    inline void core_type_convert(float &dst, int src) { dst = (float)src; }
    inline void core_type_convert(double &dst, int src) { dst = (double)src; }
    
    
    // Supported type conversions:
    
    // double               -> float
    // float                -> double
    // std::complex<double> -> std::complex<float>
    // std::complex<float>  -> std::complex<double>
    // float                -> std::complex<float>
    // double               -> std::complex<double>
    // std::complex<double> -> double (takes real part)
    // std::complex<float>  -> float (takes real part)
    // uint8_t              -> float
    // uint8_t              -> double
    // int16_t              -> float
    // int16_t              -> double
    // int                  -> float
    // int                  -> double
    
    template <typename TSRC, typename TDST> void core_block_convert(im::VecView<TDST> vdst, im::VecView<TSRC> const &vsrc);
    template <typename TSRC, typename TDST> void core_block_convert(im::MtxView<TDST> mdst, im::MtxView<TSRC> const &msrc);

}

#endif
