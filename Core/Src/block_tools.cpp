//
//  block_tools.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/25/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT> void im::core_block_clear_upper_tri(im::MtxView<TT> mdst, bool include_diagonal)
{
    IM_CHECK_VALID(mdst);
    
    int rows = mdst.rows();
    int cols = mdst.cols();
    
    int x = include_diagonal ? 0 : 1;
    
    for(int row=0; row<rows; row++)
    {
        for(int col=row+x; col<cols; col++)
            mdst(row,col) = (TT)0.0;
    }
}

#define INST(TT) template void im::core_block_clear_upper_tri(im::MtxView<TT> m, bool include_diagonal)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST
//

template <typename TT> void im::core_block_clear_lower_tri(im::MtxView<TT> mdst, bool include_diagonal)
{
    IM_CHECK_VALID(mdst);
    
    int rows = mdst.rows();
    int cols = mdst.cols();
    
    int x = include_diagonal ? 0 : 1;
    
    for(int col=0; col<cols; col++)
    {
        for(int row=col+x; row<rows; row++)
            mdst(row,col) = (TT)0.0;
    }
}

#define INST(TT) template void im::core_block_clear_lower_tri(im::MtxView<TT> m, bool include_diagonal)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST
//

template <typename TT> void im::core_block_copy_upper_tri(im::MtxView<TT> mdst, MtxView<TT> const &msrc, bool include_diagonal)
{
    IM_CHECK_VALID(msrc);
    IM_CHECK_VALID(mdst);

    int rows = std::min(msrc.rows(), msrc.cols());
    int cols = msrc.cols();
    
    IM_CHECK_ARGS(mdst.rows()>=rows && mdst.cols()>=cols);
    
    int x = include_diagonal ? 0 : 1;
    
    for(int row=0; row<rows; row++)
    {
        for(int col=row+x; col<cols; col++)
            mdst(row,col) = msrc(row,col);
    }
}

#define INST(TT) template void im::core_block_copy_upper_tri(im::MtxView<TT> mdst, MtxView<TT> const &msrc, bool include_diagonal)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST
//

template <typename TT> void im::core_block_copy_lower_tri(im::MtxView<TT> mdst, MtxView<TT> const &msrc, bool include_diagonal)
{
    IM_CHECK_VALID(msrc);
    IM_CHECK_VALID(mdst);

    int rows = msrc.rows();
    int cols = std::min(msrc.rows(), msrc.cols());
    
    IM_CHECK_ARGS(mdst.rows()>=rows && mdst.cols()>=cols);
    
    int x = include_diagonal ? 0 : 1;
    
    for(int col=0; col<cols; col++)
    {
        for(int row=col+x; row<rows; row++)
            mdst(row,col) = msrc(row,col);
    }
}

#define INST(TT) template void im::core_block_copy_lower_tri(im::MtxView<TT> mdst, MtxView<TT> const &msrc, bool include_diagonal)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_normalize(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &k)
{
    TT ss = im::core_block_blas_nrm2(vsrc);
    im::core_block_scale(vdst, vsrc, TT(1)/(ss+k));
}

#define INST(TT) template void im::core_block_normalize(VecView<TT> vdst, VecView<TT> const &vsrc, TT const &k)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_normalize(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &k)
{
    TT ss = im::core_block_blas_nrm2(msrc);
    im::core_block_scale(mdst, msrc, TT(1)/(ss+k));
}

#define INST(TT) template void im::core_block_normalize(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &k)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_normalize_columns(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &k)
{
    IM_CHECK_VALID(msrc);
    IM_CHECK_VALID(mdst);
    
    for(int i=0; i<msrc.cols(); i++)
    {
        TT ss = im::core_block_reduce_sum_squares(msrc.col(i));
        im::core_block_scale(mdst.col(i), msrc.col(i), TT(1)/(std::sqrt(ss)+k));
    }
}

#define INST(TT) template void im::core_block_normalize_columns(MtxView<TT> mdst, MtxView<TT> const &msrc, TT const &k)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//


template <typename TT> int im::core_block_nearest(MtxView<TT> const &mSet, VecView<TT> const &v)
{
    IM_CHECK_VALID(mSet);
    IM_CHECK_VALID(v);
    IM_CHECK_ARGS(mSet.rows()==v.rows());
    
    if(mSet.cols()<2)
        return 0;
    
    int nearest = 0;
    double bestdistancesq = core_real(core_block_reduce_squared_distance(mSet.col(0), v));
    
    int rows = mSet.rows();
    
    int src1stride = mSet.row_stride();
    int src2stride = v.row_stride();
    
    for(int i=1; i<mSet.cols(); i++)
    {
        TT const *psrc1 = mSet.ptr(0,i);
        TT const *psrc2 = v.ptr();
        
        int block4 = 4*(rows/4);
        int row = 0;
        
        double sumdiff = 0;
        
        for(; row<block4; row += 4, psrc1 += 4*src1stride, psrc2 += 4*src2stride)
        {
            TT d = psrc1[0] - psrc2[0];
            sumdiff += core_real(d * core_conj(d));
            
            d = psrc1[src1stride] - psrc2[src2stride];
            sumdiff += core_real(d * core_conj(d));
            
            d = psrc1[2*src1stride] - psrc2[2*src2stride];
            sumdiff += core_real(d * core_conj(d));
            
            d = psrc1[3*src1stride] - psrc2[3*src2stride];
            sumdiff += core_real(d * core_conj(d));
            
            if(sumdiff>bestdistancesq)
                goto early_exit;
        }
        
        for(; row < rows; row++, psrc1 += src1stride, psrc2 += src2stride)
        {
            TT d = psrc1[0] - psrc2[0];
            sumdiff += core_real(d * core_conj(d));
        }
        
        if(sumdiff < bestdistancesq)
        {
            bestdistancesq = sumdiff;
            nearest = i;
        }
        
    early_exit:
        ;
    }
    
    
    return nearest;
}

#define INST(TT) template int im::core_block_nearest(MtxView<TT> const &mSet, VecView<TT> const &v)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_multiply_rows(MtxView<TT> mdst, VecView<TT> const &vScale, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(msrc);
    IM_CHECK_VALID(mdst);
    IM_CHECK_VALID(vScale);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc, mdst);
    IM_CHECK_ARGS(vScale.rows()==msrc.rows());
    
    int cols = msrc.cols();
    int src_col_stride = msrc.col_stride();
    int dst_col_stride = mdst.col_stride();
    
    for(int row = 0; row<vScale.rows(); row++)
    {
        TT const *psrc = msrc.ptr(row,0);
        TT *pdst = mdst.ptr(row,0);
        TT s = vScale(row);
        
        for(int col = 0; col<cols; col++, psrc += src_col_stride, pdst += dst_col_stride)
            *pdst = *psrc * s;
    }
}

#define INST(TT) template void im::core_block_multiply_rows(MtxView<TT> mdst, VecView<TT> const &vScale, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_divide_rows(MtxView<TT> mdst, VecView<TT> const &vScale, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(msrc);
    IM_CHECK_VALID(mdst);
    IM_CHECK_VALID(vScale);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc, mdst);
    IM_CHECK_ARGS(vScale.rows()==msrc.rows());
    
    int cols = msrc.cols();
    int src_col_stride = msrc.col_stride();
    int dst_col_stride = mdst.col_stride();
    
    for(int row = 0; row<vScale.rows(); row++)
    {
        TT const *psrc = msrc.ptr(row,0);
        TT *pdst = mdst.ptr(row,0);
        TT s = (TT)1.0/vScale(row);
        
        for(int col = 0; col<cols; col++, psrc += src_col_stride, pdst += dst_col_stride)
            *pdst = *psrc * s;
    }
}

#define INST(TT) template void im::core_block_divide_rows(MtxView<TT> mdst, VecView<TT> const &vScale, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> TT im::core_block_compute_xTAx(MtxView<TT> const &mA, VecView<TT> const &vx)
{
    IM_CHECK_VALID(mA);
    IM_CHECK_VALID(vx);
    IM_CHECK_MATRIX_SQUARE(mA);
    IM_CHECK_ARGS(vx.rows()==mA.rows());
    
    int d = vx.rows();
    int strideA = mA.col_stride();
    int strideX = vx.row_stride();
    
    TT sum = (TT)0;
    for(int row=0; row<d; row++)
    {
        TT innersum = (TT)0;
        TT const *pA = mA.ptr(row,0);
        TT const *pX = vx.ptr(0);
        
        int block4 = 4*(d/4);
        int col = 0;
        for(; col<block4; col+=4, pA += 4*strideA, pX += 4*strideX)
        {
            innersum += pA[0] * pX[0];
            innersum += pA[strideA] * pX[strideX];
            innersum += pA[2*strideA] * pX[2*strideX];
            innersum += pA[3*strideA] * pX[3*strideX];
        }
        
        for(; col<d; col++, pA += strideA, pX += strideX)
            innersum += pA[0] * pX[0];
        
        sum += vx(row) * innersum;
    }
    
    return sum;
}

#define INST(TT) template TT im::core_block_compute_xTAx(MtxView<TT> const &mA, VecView<TT> const &vx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST



