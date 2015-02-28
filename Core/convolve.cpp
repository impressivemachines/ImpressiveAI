//
//  convolve.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/3/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::MatrixPadder<TT>::pad_zero(MtxRect const &rct)
{
    int const col_stride = m_block_store.col_stride();
    
    for(int r=0; r<rct.size.rows; r++)
    {
        int sr = r + rct.origin.row;
        TT *p = m_block_store.ptr(r,0);
        
        if(sr<0 || sr>=m_source.rows())
        {
            // top/bottom pad
            for(int c=0; c<rct.size.cols; c++, p += col_stride)
                *p = (TT)0;
        }
        else
        {
            for(int c=0; c<rct.size.cols; c++, p += col_stride)
            {
                int sc = c + rct.origin.col;
                
                if(sc<0 || sc>=m_source.cols())
                    *p = (TT)0;
                else
                    *p = m_source(sr,sc);
            }
        }
    }
}

template <typename TT>
void im::MatrixPadder<TT>::pad_extend(MtxRect const &rct)
{
    int const col_stride = m_block_store.col_stride();
    
    for(int r=0; r<rct.size.rows; r++)
    {
        int sr = r + rct.origin.row;
        sr = core_clamp(sr, m_source.rows());
        
        TT *p = m_block_store.ptr(r,0);
        for(int c=0; c<rct.size.cols; c++, p += col_stride)
        {
            int sc = c + rct.origin.col;
            sc = core_clamp(sc, m_source.cols());
            
            *p = m_source(sr,sc);
        }
    }
}

template <typename TT>
void im::MatrixPadder<TT>::pad_wrap(MtxRect const &rct)
{
    int const col_stride = m_block_store.col_stride();
    
    for(int r=0; r<rct.size.rows; r++)
    {
        int sr = r + rct.origin.row;
        sr = core_modulus(sr, m_source.rows());
        
        TT *p = m_block_store.ptr(r,0);
        for(int c=0; c<rct.size.cols; c++, p += col_stride)
        {
            int sc = c + rct.origin.col;
            sc = core_modulus(sc, m_source.cols());
            
            *p = m_source(sr,sc);
        }
    }
}

template <typename TT>
void im::MatrixPadder<TT>::pad_reflect(MtxRect const &rct)
{
    int const col_stride = m_block_store.col_stride();
    
    for(int r=0; r<rct.size.rows; r++)
    {
        int sr = r + rct.origin.row;
        sr = core_reflect(sr, m_source.rows());
        
        TT *p = m_block_store.ptr(r,0);
        for(int c=0; c<rct.size.cols; c++, p += col_stride)
        {
            int sc = c + rct.origin.col;
            sc = core_reflect(sc, m_source.cols());
            
            *p = m_source(sr,sc);
        }
    }
}

template <typename TT>
im::MtxView<TT> im::MatrixPadder<TT>::get(MtxRect const &rct, PadMode pad)
{
    int rows = rct.size.rows;
    int cols = rct.size.cols;
    
    MtxRect rctsrc(0,0,m_source.rows(),m_source.cols());
    
    if(rct.is_within(rctsrc))
        return m_source.block(rct);
    
    if(m_block_store.rows() < rows || m_block_store.cols() < cols)
        m_block_store.resize(rows,cols);
    
    switch(pad)
    {
        case PadModeZero:
        default:
            pad_zero(rct);
        break;
      
        case PadModeExtend:
            pad_extend(rct);
        break;
            
        case PadModeWrap:
            pad_wrap(rct);
        break;
            
        case PadModeReflect:
            pad_reflect(rct);
        break;
    }
    
    return m_block_store.view().block(0,0,rows,cols);
}

template class im::MatrixPadder<float>;
template class im::MatrixPadder<double>;
template class im::MatrixPadder<im::Cf>;
template class im::MatrixPadder<im::Cd>;


// General 2D convolution
// The kernel center indicates the zero origin of the kernel matrix
// The start location indicates the position in the input matrix corresponding to the top left of the output matrix
template <typename TT>
void im::core_convolve(MtxView<TT> mavout, MtxView<TT> const &mavkernel, MtxView<TT> const &mavin,
                          PadMode pad, MtxLoc kernel_center, MtxLoc start_loc, TT output_pre_weight, int drows, int dcols, int blocksize)
{
    IM_CHECK_VALID(mavin);
    IM_CHECK_VALID(mavout);
    IM_CHECK_VALID(mavkernel);
    IM_CHECK_ARGS(kernel_center.row>=0 && kernel_center.col>=0);
    IM_CHECK_ARGS(kernel_center.row<mavkernel.rows() && kernel_center.col<mavkernel.cols());
    IM_CHECK_ARGS(drows>0);
    IM_CHECK_ARGS(dcols>0);
    IM_CHECK_ARGS(blocksize>0);
    
    // Switch to correlation
    MtxView<TT> mkern = mavkernel.reverse_cols().reverse_rows();
    kernel_center.row = mavkernel.rows() - 1 - kernel_center.row;
    kernel_center.col = mavkernel.cols() - 1 - kernel_center.col;
    
    int krows = mkern.rows();
    int kcols = mkern.cols();
    
    RectSubdivider rsd;
    MatrixPadder<TT> mp;

    rsd.init(MtxRect(0,0,mavout.rows(),mavout.cols()), MtxSize(blocksize, blocksize));
    mp.init(mavin);

    MtxRect dstrct;
    while(rsd.next(dstrct))
    {
        MtxRect srcrct;
        srcrct.origin.row = dstrct.origin.row * drows - kernel_center.row + start_loc.row;
        srcrct.origin.col = dstrct.origin.col * dcols - kernel_center.col + start_loc.col;
        srcrct.size.rows = (dstrct.size.rows - 1) * drows + krows;
        srcrct.size.cols = (dstrct.size.cols - 1) * dcols + kcols;
        
       // printf("dst %d %d %d %d\n", dstrct.origin.row, dstrct.origin.col, dstrct.size.rows, dstrct.size.cols);
       // printf("src %d %d %d %d\n", srcrct.origin.row, srcrct.origin.col, srcrct.size.rows, srcrct.size.cols);
        
        MtxView<TT> msrc = mp.get(srcrct, pad);
        MtxView<TT> mdst = mavout.block(dstrct);
       // msrc.print();
        
        if(output_pre_weight==(TT)0)
            mdst = (TT)0;
        else if(output_pre_weight!=(TT)1)
            core_block_scale(mdst, mdst, output_pre_weight);
        
        if(kcols>=krows)
        {
            for(int i=0; i<dstrct.size.rows; i++)
                for(int j=0; j<dstrct.size.cols; j++)
                {
                    TT const *psrc1 = msrc.ptr(i*drows, j*dcols);
                    TT const *pk1 = mkern.ptr();
                    
                    TT sum = (TT)0;
                    
                    for(int kr = 0; kr<krows; kr++, psrc1 += msrc.row_stride(), pk1 += mkern.row_stride())
                        core_madd_loop(sum, psrc1, msrc.col_stride(), pk1, mkern.col_stride(), kcols);
                    
                    mdst(i, j) += sum;
                }
        }
        else
        {
            for(int i=0; i<dstrct.size.rows; i++)
                for(int j=0; j<dstrct.size.cols; j++)
                {
                    TT const *psrc1 = msrc.ptr(i*drows, j*dcols);
                    TT const *pk1 = mkern.ptr();
                    
                    TT sum = (TT)0;
                    
                    for(int kc = 0; kc<kcols; kc++, psrc1 += msrc.col_stride(), pk1 += mkern.col_stride())
                        core_madd_loop(sum, psrc1, msrc.row_stride(), pk1, mkern.row_stride(), krows);
                    
                    mdst(i, j) += sum;
                }
        }
    }
}

#define INST(TT) template void im::core_convolve(MtxView<TT> mavout, MtxView<TT> const &mavkernel, MtxView<TT> const &mavin, PadMode pad, MtxLoc kernel_center, MtxLoc start_loc, TT output_pre_weight, int drows, int dcols, int blocksize)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

// Specialization for 1D
// The kernel center indicates the zero origin of the kernel matrix
// The start location indicates the position in the input matrix corresponding to the first sample of the output matrix
template <typename TT>
void im::core_convolve(VecView<TT> vvout, VecView<TT> const &vvkernel, VecView<TT> const &vvin,
                          PadMode pad, int kernel_center, int start_loc, TT output_pre_weight, int drows)
{
    IM_CHECK_VALID(vvin);
    IM_CHECK_VALID(vvout);
    IM_CHECK_VALID(vvkernel);
    IM_CHECK_ARGS(kernel_center>=0 && kernel_center<vvkernel.rows());
    IM_CHECK_ARGS(drows>0);
    
    int size_k = vvkernel.rows();
    int size_in = vvin.rows();
    int size_out = vvout.rows();
    
    // switch to correlation
    VecView<TT> vvk = vvkernel.reverse();
    kernel_center = size_k - 1 - kernel_center;
    
    if(output_pre_weight==(TT)0)
        vvout = (TT)0;
    else if(output_pre_weight!=(TT)1)
        core_block_scale(vvout, vvout, output_pre_weight);
    
    int ks = start_loc - kernel_center;
    int upper = size_in - size_k;
    
    int kstride = vvk.row_stride();
    int instride = vvin.row_stride();
    
    for(int i=0; i<size_out; i++, ks += drows)
    {
        TT sum = (TT)0;

        if(ks >= 0 && ks <= upper)
        {
            // no padding
            core_madd_loop(sum, vvk.ptr(), kstride, vvin.ptr(ks), instride, size_k);
        }
        else
        {
            // padding
            switch(pad)
            {
                case im::PadModeZero:
                default:
                    for(int k=0; k<size_k; k++)
                    {
                        int rink = ks + k;
                        if(rink >= 0 && rink < size_in)
                            sum += vvk(k) * vvin(rink);
                    }
                    break;
                    
                case im::PadModeExtend:
                    for(int k=0; k<size_k; k++)
                        sum += vvk(k) * vvin(im::core_clamp(ks + k, size_in));
                    break;
                    
                case im::PadModeWrap:
                    for(int k=0; k<size_k; k++)
                        sum += vvk(k) * vvin(im::core_modulus(ks + k, size_in));
                    break;
                    
                case im::PadModeReflect:
                    for(int k=0; k<size_k; k++)
                        sum += vvk(k) * vvin(im::core_reflect(ks + k, size_in));
                    break;
            }
        }
        
        vvout(i) += sum;
    }
}

#define INST(TT) template void im::core_convolve(VecView<TT> vvout, VecView<TT> const &vvkernel, VecView<TT> const &vvin, PadMode pad, int kernel_center, int start_loc, TT output_pre_weight, int drows)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

