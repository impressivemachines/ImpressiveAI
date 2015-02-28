//
//  block_convert.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 1/25/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TSRC, typename TDST> void im::core_block_convert(im::VecView<TDST> vdst, im::VecView<TSRC> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    IM_CHECK_VALID(vdst);
    IM_CHECK_VECTOR_SIZES_MATCH(vsrc,vdst);
    
    int rows = vsrc.rows();
    int srcstride = vsrc.row_stride();
    int dststride = vdst.row_stride();
    
    TSRC const *psrc = vsrc.ptr();
    TDST *pdst = vdst.ptr();
    
    int block4 = 4*(rows/4);
    int row = 0;
    for(; row < block4; row += 4, pdst += 4*dststride, psrc += 4*srcstride)
    {
        core_type_convert(pdst[0], psrc[0]);
        core_type_convert(pdst[dststride], psrc[srcstride]);
        core_type_convert(pdst[2*dststride], psrc[2*srcstride]);
        core_type_convert(pdst[3*dststride], psrc[3*srcstride]);
    }
    
    for(; row < rows; row++, pdst += dststride, psrc += srcstride)
        core_type_convert(pdst[0], psrc[0]);
}

template void im::core_block_convert(im::VecView<double> vdst, im::VecView<float> const &vsrc);
template void im::core_block_convert(im::VecView<float> vdst, im::VecView<double> const &vsrc);
template void im::core_block_convert(im::VecView<im::Cd> vdst, im::VecView<im::Cf> const &vsrc);
template void im::core_block_convert(im::VecView<im::Cf> vdst, im::VecView<im::Cd> const &vsrc);
template void im::core_block_convert(im::VecView<im::Cf> vdst, im::VecView<float> const &vsrc);
template void im::core_block_convert(im::VecView<im::Cd> vdst, im::VecView<double> const &vsrc);
template void im::core_block_convert(im::VecView<float> vdst, im::VecView<im::Cf> const &vsrc);
template void im::core_block_convert(im::VecView<double> vdst, im::VecView<im::Cd> const &vsrc);

template void im::core_block_convert(im::VecView<float> vdst, im::VecView<uint8_t> const &vsrc);
template void im::core_block_convert(im::VecView<float> vdst, im::VecView<int16_t> const &vsrc);
template void im::core_block_convert(im::VecView<float> vdst, im::VecView<int> const &vsrc);
template void im::core_block_convert(im::VecView<double> vdst, im::VecView<uint8_t> const &vsrc);
template void im::core_block_convert(im::VecView<double> vdst, im::VecView<int16_t> const &vsrc);
template void im::core_block_convert(im::VecView<double> vdst, im::VecView<int> const &vsrc);

template <typename TSRC, typename TDST> void im::core_block_convert(im::MtxView<TDST> mdst, im::MtxView<TSRC> const &msrc)
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
            TSRC const *psrc = msrc.ptr(0,col);
            TDST *pdst = mdst.ptr(0,col);
            
            int block4 = 4*(rows/4);
            int row = 0;
            for(; row < block4; row += 4, pdst += 4*dststride, psrc += 4*srcstride)
            {
                core_type_convert(pdst[0], psrc[0]);
                core_type_convert(pdst[dststride], psrc[srcstride]);
                core_type_convert(pdst[2*dststride], psrc[2*srcstride]);
                core_type_convert(pdst[3*dststride], psrc[3*srcstride]);
            }
            
            for(; row < rows; row++, pdst += dststride, psrc += srcstride)
                core_type_convert(pdst[0], psrc[0]);
        }
    }
    else
    {
        int dststride = mdst.col_stride(); // note that strides can be negative
        int srcstride = msrc.col_stride();
        
        for(int row = 0; row < rows; row++)
        {
            TSRC const *psrc = msrc.ptr(row,0);
            TDST *pdst = mdst.ptr(row,0);
            
            int block4 = 4*(cols/4);
            int col = 0;
            for(; col < block4; col += 4, pdst += 4*dststride, psrc += 4*srcstride)
            {
                core_type_convert(pdst[0], psrc[0]);
                core_type_convert(pdst[dststride], psrc[srcstride]);
                core_type_convert(pdst[2*dststride], psrc[2*srcstride]);
                core_type_convert(pdst[3*dststride], psrc[3*srcstride]);
            }
            
            for(; col < cols; col++, pdst += dststride, psrc += srcstride)
                core_type_convert(pdst[0], psrc[0]);
        }
    }
}

template void im::core_block_convert(im::MtxView<double> mdst, im::MtxView<float> const &msrc);
template void im::core_block_convert(im::MtxView<float> mdst, im::MtxView<double> const &msrc);
template void im::core_block_convert(im::MtxView<im::Cd> mdst, im::MtxView<im::Cf> const &msrc);
template void im::core_block_convert(im::MtxView<im::Cf> mdst, im::MtxView<im::Cd> const &msrc);
template void im::core_block_convert(im::MtxView<im::Cf> mdst, im::MtxView<float> const &msrc);
template void im::core_block_convert(im::MtxView<im::Cd> mdst, im::MtxView<double> const &msrc);
template void im::core_block_convert(im::MtxView<float> mdst, im::MtxView<im::Cf> const &msrc);
template void im::core_block_convert(im::MtxView<double> mdst, im::MtxView<im::Cd> const &msrc);

template void im::core_block_convert(im::MtxView<float> mdst, im::MtxView<uint8_t> const &msrc);
template void im::core_block_convert(im::MtxView<float> mdst, im::MtxView<int16_t> const &msrc);
template void im::core_block_convert(im::MtxView<float> mdst, im::MtxView<int> const &msrc);
template void im::core_block_convert(im::MtxView<double> mdst, im::MtxView<uint8_t> const &msrc);
template void im::core_block_convert(im::MtxView<double> mdst, im::MtxView<int16_t> const &msrc);
template void im::core_block_convert(im::MtxView<double> mdst, im::MtxView<int> const &msrc);


