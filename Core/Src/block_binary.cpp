//
//  block_binary.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT> struct BlendParam
{
    TT ax;
    TT max;
};

template <typename TT, typename P>
struct BinaryAdd
{
    static inline TT apply(TT const &x, TT const &y, P const &param) { return x + y; }
};

template <typename TT, typename P>
struct BinarySub
{
    static inline TT apply(TT const &x, TT const &y, P const &param) { return x - y; }
};

template <typename TT, typename P>
struct BinaryMultiply
{
    static inline TT apply(TT const &x, TT const &y, P const &param) { return x * y; }
};

template <typename TT, typename P>
struct BinaryDivide
{
    static inline TT apply(TT const &x, TT const &y, P const &param) { return x / y; }
};

template <typename TT, typename P>
struct BinaryBlend
{
    static inline TT apply(TT const &x, TT const &y, P const &param) { return x * param.max + y * param.ax; }
};

template <typename TT, typename P>
struct BinaryAddScaled
{
    static inline TT apply(TT const &x, TT const &y, P const &param) { return x + y * param; }
};

template <typename TT, typename OP, typename P> void core_priv_block_binary_op(im::VecView<TT> vdst, im::VecView<TT> const &vsrc1, im::VecView<TT> const &vsrc2, P const &param)
{
    IM_CHECK_VALID(vsrc1);
    IM_CHECK_VALID(vsrc2);
    IM_CHECK_VALID(vdst);
    IM_CHECK_VECTOR_SIZES_MATCH(vsrc1,vsrc2);
    IM_CHECK_VECTOR_SIZES_MATCH(vsrc1,vdst);
    
    // TODO: add SSE here
    
    // Fallback implementation
    int rows = vsrc1.rows();
    int block4 = 4*(rows/4);
    int src1stride = vsrc1.row_stride(); // note that strides can be negative
    int src2stride = vsrc2.row_stride();
    int dststride = vdst.row_stride();
    
    TT const *psrc1 = vsrc1.ptr();
    TT const *psrc2 = vsrc2.ptr();
    TT *pdst = vdst.ptr();
    
    int row = 0;
    
    for(; row < block4; row+=4, pdst += 4*dststride, psrc1 += 4*src1stride, psrc2 += 4*src2stride)
    {
        pdst[0] = OP::apply(psrc1[0], psrc2[0], param);
        pdst[dststride] = OP::apply(psrc1[src1stride], psrc2[src2stride], param);
        pdst[2*dststride] = OP::apply(psrc1[2*src1stride], psrc2[2*src2stride], param);
        pdst[3*dststride] = OP::apply(psrc1[3*src1stride], psrc2[3*src2stride], param);
    }
    
    for(; row < rows; row++, pdst += dststride, psrc1 += src1stride, psrc2 += src2stride)
        pdst[0] = OP::apply(psrc1[0], psrc2[0], param);
}

template <typename TT> void im::core_block_add_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    core_priv_block_binary_op<TT, BinaryAdd<TT, float>, float>(vdst, vsrc1, vsrc2, 0.0f);
}

#define INST(TT) template void im::core_block_add_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sub_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    core_priv_block_binary_op<TT, BinarySub<TT, float>, float>(vdst, vsrc1, vsrc2, 0.0f);
}

#define INST(TT) template void im::core_block_sub_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_multiply_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    core_priv_block_binary_op<TT, BinaryMultiply<TT, float>, float>(vdst, vsrc1, vsrc2, 0.0f);
}

#define INST(TT) template void im::core_block_multiply_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_divide_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    core_priv_block_binary_op<TT, BinaryDivide<TT, float>, float>(vdst, vsrc1, vsrc2, 0.0f);
}

#define INST(TT) template void im::core_block_divide_elements(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_blend(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2, float alpha)
{
    BlendParam<TT> bp;
    bp.ax = (TT)alpha;
    bp.max = (TT)(1-alpha);
    core_priv_block_binary_op<TT, BinaryBlend<TT, BlendParam<TT> >, BlendParam<TT> >(vdst, vsrc1, vsrc2, bp);
}

#define INST(TT) template void im::core_block_blend(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2, float alpha)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_add_scaled(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2, TT const &scale)
{
    core_priv_block_binary_op<TT, BinaryAddScaled<TT, TT>, TT>(vdst, vsrc1, vsrc2, scale);
}

#define INST(TT) template void im::core_block_add_scaled(VecView<TT> vdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2, TT const &scale)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST


template <typename TT, typename OP, typename P> void core_priv_block_binary_op(im::MtxView<TT> mdst, im::MtxView<TT> const &msrc1, im::MtxView<TT> const &msrc2, P const &param)
{
    IM_CHECK_VALID(msrc1);
    IM_CHECK_VALID(msrc2);
    IM_CHECK_VALID(mdst);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc1,msrc2);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc1,mdst);
    
    // TODO: add SSE here
    
    // Fallback implementation
    int rows = msrc1.rows();
    int cols = msrc1.cols();
    
    if(rows>cols)
    {
        int src1stride = msrc1.row_stride(); // note that strides can be negative
        int src2stride = msrc2.row_stride();
        int dststride = mdst.row_stride();
        
        for(int col = 0; col < cols; col++)
        {
            TT const *psrc1 = msrc1.ptr(0,col);
            TT const *psrc2 = msrc2.ptr(0,col);
            TT *pdst = mdst.ptr(0,col);
            
            int block4 = 4*(rows/4);
            int row = 0;
            for(; row < block4; row += 4, pdst += 4*dststride, psrc1 += 4*src1stride, psrc2 += 4*src2stride)
            {
                pdst[0] = OP::apply(psrc1[0], psrc2[0], param);
                pdst[dststride] = OP::apply(psrc1[src1stride], psrc2[src2stride], param);
                pdst[2*dststride] = OP::apply(psrc1[2*src1stride], psrc2[2*src2stride], param);
                pdst[3*dststride] = OP::apply(psrc1[3*src1stride], psrc2[3*src2stride], param);
            }
            
            for(; row < rows; row++, pdst += dststride, psrc1 += src1stride, psrc2 += src2stride)
                pdst[0] = OP::apply(psrc1[0], psrc2[0], param);
        }
    }
    else
    {
        int src1stride = msrc1.col_stride(); // note that strides can be negative
        int src2stride = msrc2.col_stride();
        int dststride = mdst.col_stride();
        
        for(int row = 0; row < rows; row++)
        {
            TT const *psrc1 = msrc1.ptr(row,0);
            TT const *psrc2 = msrc2.ptr(row,0);
            TT *pdst = mdst.ptr(row,0);
            
            int block4 = 4*(cols/4);
            int col = 0;
            
            for(; col<block4; col += 4, pdst += 4*dststride, psrc1 += 4*src1stride, psrc2 += 4*src2stride)
            {
                pdst[0] = OP::apply(psrc1[0], psrc2[0], param);
                pdst[dststride] = OP::apply(psrc1[src1stride], psrc2[src2stride], param);
                pdst[2*dststride] = OP::apply(psrc1[2*src1stride], psrc2[2*src2stride], param);
                pdst[3*dststride] = OP::apply(psrc1[3*src1stride], psrc2[3*src2stride], param);
            }
            
            for(; col < cols; col++, pdst += dststride, psrc1 += src1stride, psrc2 += src2stride)
                pdst[0] = OP::apply(psrc1[0], psrc2[0], param);
        }
    }
}

template <typename TT> void im::core_block_add_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
{
    core_priv_block_binary_op<TT, BinaryAdd<TT, float>, float>(mdst, msrc1, msrc2, 0.0f);
}

#define INST(TT) template void im::core_block_add_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sub_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
{
    core_priv_block_binary_op<TT, BinarySub<TT, float>, float>(mdst, msrc1, msrc2, 0.0f);
}

#define INST(TT) template void im::core_block_sub_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_multiply_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
{
    core_priv_block_binary_op<TT, BinaryMultiply<TT, float>, float>(mdst, msrc1, msrc2, 0.0f);
}

#define INST(TT) template void im::core_block_multiply_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_divide_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
{
    core_priv_block_binary_op<TT, BinaryDivide<TT, float>, float>(mdst, msrc1, msrc2, 0.0f);
}

#define INST(TT) template void im::core_block_divide_elements(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_blend(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2, float alpha)
{
    BlendParam<TT> bp;
    bp.ax = (TT)alpha;
    bp.max = (TT)(1-alpha);
    core_priv_block_binary_op<TT, BinaryBlend<TT, BlendParam<TT> >, BlendParam<TT> >(mdst, msrc1, msrc2, bp);
}

#define INST(TT) template void im::core_block_blend(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2, float alpha)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_add_scaled(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2, TT const &scale)
{
    core_priv_block_binary_op<TT, BinaryAddScaled<TT, TT>, TT>(mdst, msrc1, msrc2, scale);
}

#define INST(TT) template void im::core_block_add_scaled(MtxView<TT> mdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2, TT const &scale)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST


template <typename TT> void im::core_block_outer_product(MtxView<TT> mdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    IM_CHECK_VALID(vsrc1);
    IM_CHECK_VALID(vsrc2);
    IM_CHECK_VALID(mdst);
    IM_CHECK_ARGS(mdst.rows()==vsrc1.rows());
    IM_CHECK_ARGS(mdst.cols()==vsrc2.rows());
    
    int rows = mdst.rows();
    int cols = mdst.cols();
    int mdststride = mdst.col_stride();
    int src2stride = vsrc2.row_stride();
    
    for(int row = 0; row < rows; row++)
    {
        TT v = vsrc1(row);
        TT const *psrc2 = vsrc2.ptr(0);
        TT *pdst = mdst.ptr(row,0);
        int col = 0;
        int block4 = 4*(cols/4);
        
        for(; col < block4; col+=4, psrc2 += 4*src2stride, pdst += 4*mdststride)
        {
            pdst[0] = v * psrc2[0];
            pdst[mdststride] = v * psrc2[src2stride];
            pdst[2*mdststride] = v * psrc2[2*src2stride];
            pdst[3*mdststride] = v * psrc2[3*src2stride];
        }
        
        for(; col < cols; col++, psrc2 += src2stride, pdst += mdststride)
            pdst[0] = v * psrc2[0];
    }
}

#define INST(TT) template void im::core_block_outer_product(MtxView<TT> mdst, VecView<TT> const &vsrc1, VecView<TT> const &vsrc2);
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST


