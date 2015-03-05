//
//  block_unary.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT> struct ScaleOffsetParam
{
    TT scale;
    TT offset;
};

//template <typename TT, typename P>
//struct UnaryCopy
//{
//    static inline TT apply(TT const &x, P const &param) { return x; }
//};

template <typename TT, typename P>
struct UnaryNeg
{
    static inline TT apply(TT const &x, P const &param) { return -x; }
};

template <typename TT, typename P>
struct UnaryConj
{
    static inline TT apply(TT const &x, P const &param) { return im::core_conj(x); }
};

template <typename TT, typename P>
struct UnaryRecip
{
    static inline TT apply(TT const &x, P const &param) { return (TT)1 / (x + param); }
};

template <typename TT, typename P>
struct UnaryScale
{
    static inline TT apply(TT const &x, P const &param) { return x * param; }
};

template <typename TT, typename P>
struct UnaryOffset
{
    static inline TT apply(TT const &x, P const &param) { return x + param; }
};

template <typename TT, typename P>
struct UnaryScaleOffset
{
    static inline TT apply(TT const &x, P const &param) { return x * param.scale + param.offset; }
};

template <typename TT, typename P>
struct UnarySqrt
{
    static inline TT apply(TT const &x, P const &param) { return std::sqrt(x); }
};

template <typename TT, typename P>
struct UnaryAbs
{
    static inline TT apply(TT const &x, P const &param) { return std::abs(x); }
};

template <typename TT, typename P>
struct UnaryExp
{
    static inline TT apply(TT const &x, P const &param) { return std::exp(x); }
};

template <typename TT, typename P>
struct UnaryLog
{
    static inline TT apply(TT const &x, P const &param) { return std::log(x); }
};

template <typename TT, typename P>
struct UnaryPow
{
    static inline TT apply(TT const &x, P const &param) { return std::pow(x, param); }
};

template <typename TT, typename P>
struct UnarySigm
{
    static inline TT apply(TT const &x, P const &param) { return (TT)1/((TT)1 + std::exp(-x)); }
};


template <typename TT, typename OP, typename P> void core_priv_block_unary_op(im::VecView<TT> vdst, im::VecView<TT> const &vsrc, P const &param)
{
    IM_CHECK_VALID(vsrc);
    IM_CHECK_VALID(vdst);
    IM_CHECK_VECTOR_SIZES_MATCH(vsrc,vdst);
    
    // TODO: add SSE here
    
    // Fallback implementation
    int rows = vdst.rows();
    int block4 = 4*(rows/4);
    int dststride = vdst.row_stride(); // note that strides can be negative
    int srcstride = vsrc.row_stride();
    
    TT const *psrc = vsrc.ptr();
    TT *pdst = vdst.ptr();
    
    int row = 0;
    
    for(; row < block4; row+=4, pdst += 4*dststride, psrc += 4*srcstride)
    {
        pdst[0] = OP::apply(psrc[0], param);
        pdst[dststride] = OP::apply(psrc[srcstride], param);
        pdst[2*dststride] = OP::apply(psrc[2*srcstride], param);
        pdst[3*dststride] = OP::apply(psrc[3*srcstride], param);
    }
    
    for(; row < rows; row++, pdst += dststride, psrc += srcstride)
        pdst[0] = OP::apply(psrc[0], param);
}

template <typename TT> void im::core_block_neg(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnaryNeg<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_neg(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_conj(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnaryConj<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_conj(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_recip(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &epsilon)
{
    core_priv_block_unary_op<TT, UnaryRecip<TT, TT>, TT>(vdst, vsrc, epsilon);
}

#define INST(TT) template void im::core_block_recip(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &epsilon)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_scale(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale)
{
    core_priv_block_unary_op<TT, UnaryScale<TT, TT>, TT>(vdst, vsrc, scale);
}

#define INST(TT) template void im::core_block_scale(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_offset(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &offset)
{
    core_priv_block_unary_op<TT, UnaryOffset<TT, TT>, TT>(vdst, vsrc, offset);
}

#define INST(TT) template void im::core_block_offset(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_scale_offset(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale, TT const &offset)
{
    ScaleOffsetParam<TT> so;
    so.scale = scale;
    so.offset = offset;
    
    core_priv_block_unary_op<TT, UnaryScaleOffset<TT, ScaleOffsetParam<TT>>, ScaleOffsetParam<TT>>(vdst, vsrc, so);
}

#define INST(TT) template void im::core_block_scale_offset(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &scale, TT const &offset)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_abs(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnaryAbs<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_abs(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sqrt(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnarySqrt<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_sqrt(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_pow(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &exponent)
{
    core_priv_block_unary_op<TT, UnaryPow<TT, TT>, TT>(vdst, vsrc, exponent);
}

#define INST(TT) template void im::core_block_pow(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &exponent)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_exp(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnaryExp<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_exp(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_log(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnaryLog<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_log(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sigm(VecView<TT> vdst, VecView<TT> const &vsrc)
{
    core_priv_block_unary_op<TT, UnarySigm<TT, float>, float>(vdst, vsrc, 0.0f);
}

#define INST(TT) template void im::core_block_sigm(VecView<TT> vdst, VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT, typename OP, typename P> void core_priv_block_unary_op(im::MtxView<TT> mdst, im::MtxView<TT> const &msrc, P const &param)
{
    IM_CHECK_VALID(msrc);
    IM_CHECK_VALID(mdst);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc,mdst);
    
    // TODO: add SSE here
    
    // Fallback implementation
    int rows = mdst.rows();
    int cols = mdst.cols();
    
    if(rows>cols)
    {
        int dststride = mdst.row_stride(); // note that strides can be negative
        int srcstride = msrc.row_stride();
        
        for(int col = 0; col < cols; col++)
        {
            TT const *psrc = msrc.ptr(0,col);
            TT *pdst = mdst.ptr(0,col);
            
            int block4 = 4*(rows/4);
            int row = 0;
            for(; row < block4; row += 4, pdst += 4*dststride, psrc += 4*srcstride)
            {
                pdst[0] = OP::apply(psrc[0], param);
                pdst[dststride] = OP::apply(psrc[srcstride], param);
                pdst[2*dststride] = OP::apply(psrc[2*srcstride], param);
                pdst[3*dststride] = OP::apply(psrc[3*srcstride], param);
            }
            
            for(; row < rows; row++, pdst += dststride, psrc += srcstride)
                pdst[0] = OP::apply(psrc[0], param);
        }
    }
    else
    {
        int dststride = mdst.col_stride(); // note that strides can be negative
        int srcstride = msrc.col_stride();
        
        for(int row = 0; row < rows; row++)
        {
            TT const *psrc = msrc.ptr(row,0);
            TT *pdst = mdst.ptr(row,0);
            
            int block4 = 4*(cols/4);
            int col = 0;
            for(; col < block4; col += 4, pdst += 4*dststride, psrc += 4*srcstride)
            {
                pdst[0] = OP::apply(psrc[0], param);
                pdst[dststride] = OP::apply(psrc[srcstride], param);
                pdst[2*dststride] = OP::apply(psrc[2*srcstride], param);
                pdst[3*dststride] = OP::apply(psrc[3*srcstride], param);
            }
            
            for(; col < cols; col++, pdst += dststride, psrc += srcstride)
                pdst[0] = OP::apply(psrc[0], param);
        }
    }
}

template <typename TT> void im::core_block_neg(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnaryNeg<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_neg(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_conj(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnaryConj<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_conj(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_recip(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &epsilon)
{
    core_priv_block_unary_op<TT, UnaryRecip<TT, TT>, TT>(mdst, msrc, epsilon);
}

#define INST(TT) template void im::core_block_recip(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &epsilon)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_scale(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale)
{
    core_priv_block_unary_op<TT, UnaryScale<TT, TT>, TT>(mdst, msrc, scale);
}

#define INST(TT) template void im::core_block_scale(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_offset(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &offset)
{
    core_priv_block_unary_op<TT, UnaryOffset<TT, TT>, TT>(mdst, msrc, offset);
}

#define INST(TT) template void im::core_block_offset(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_scale_offset(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale, TT const &offset)
{
    ScaleOffsetParam<TT> so;
    so.scale = scale;
    so.offset = offset;
    
    core_priv_block_unary_op<TT, UnaryScaleOffset<TT, ScaleOffsetParam<TT>>, ScaleOffsetParam<TT>>(mdst, msrc, so);
}

#define INST(TT) template void im::core_block_scale_offset(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &scale, TT const &offset)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_abs(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnaryAbs<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_abs(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sqrt(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnarySqrt<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_sqrt(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_pow(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &exponent)
{
    core_priv_block_unary_op<TT, UnaryPow<TT, TT>, TT>(mdst, msrc, exponent);
}

#define INST(TT) template void im::core_block_pow(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &exponent)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_exp(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnaryExp<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_exp(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_log(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnaryLog<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_log(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sigm(MtxView<TT> mdst, MtxView<TT> const &msrc)
{
    core_priv_block_unary_op<TT, UnarySigm<TT, float>, float>(mdst, msrc, 0.0f);
}

#define INST(TT) template void im::core_block_sigm(MtxView<TT> mdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST




