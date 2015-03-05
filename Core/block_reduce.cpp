//
//  block_reduce.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
struct ReduceAdd
{
    static inline TT apply(TT const &x, TT const &tally) { return x + tally; }
};

template <typename TT>
struct ReduceMultiply
{
    static inline TT apply(TT const &x, TT const &tally) { return x * tally; }
};

template <typename TT>
struct ReduceMax
{
    static inline TT apply(TT const &x, TT const &tally) { return (TT)std::max(im::core_real(x), im::core_real(tally)); }
};

template <typename TT>
struct ReduceMin
{
    static inline TT apply(TT const &x, TT const &tally) { return (TT)std::min(im::core_real(x), im::core_real(tally)); }
};

template <typename TT>
struct ReduceMaxAbs
{
    static inline TT apply(TT const &x, TT const &tally) { return (TT)std::max(std::abs(x), std::abs(tally)); }
};

template <typename TT>
struct ReduceMinAbs
{
    static inline TT apply(TT const &x, TT const &tally) { return (TT)std::min(std::abs(x), std::abs(tally)); }
};

template <typename TT>
struct ReduceSumSquares
{
    static inline TT apply(TT const &x, TT const &tally) { return x * im::core_conj(x) + tally; }
};

template <typename TT, typename OP> TT core_priv_block_reduce_op(im::VecView<TT> const &vsrc, TT const &tally_start)
{
    IM_CHECK_VALID(vsrc);
    
    TT tally = tally_start;
    
    int rows = vsrc.rows();
    int srcstride = vsrc.row_stride(); // note that strides can be negative
    int block4 = 4*(rows/4);
    int row = 0;
    TT const *psrc = vsrc.ptr();
    
    for(; row < block4; row+=4, psrc += 4*srcstride)
    {
        tally = OP::apply(psrc[0], tally);
        tally = OP::apply(psrc[srcstride], tally);
        tally = OP::apply(psrc[2*srcstride], tally);
        tally = OP::apply(psrc[3*srcstride], tally);
    }
    
    for(; row < rows; row++, psrc += srcstride)
        tally = OP::apply(psrc[0], tally);

    return tally;
}

template <typename TT> TT im::core_block_reduce_add(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceAdd<TT>>(vsrc, (TT)0);
}

#define INST(TT) template TT im::core_block_reduce_add(VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_multiply(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceMultiply<TT>>(vsrc, (TT)1);
}

#define INST(TT) template TT im::core_block_reduce_multiply(VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_mean(VecView<TT> const &vsrc)
{
    if(vsrc.count()<1)
        return (TT)0;
    
    return core_priv_block_reduce_op<TT, ReduceAdd<TT>>(vsrc, (TT)0) / (TT)vsrc.count();
}

#define INST(TT) template TT im::core_block_reduce_mean(VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_variance(VecView<TT> const &vsrc, TT *pmeanrtn)
{
    IM_CHECK_VALID(vsrc);
    
    if(vsrc.count()<1)
    {
        if(pmeanrtn)
            *pmeanrtn = (TT)0;
        return (TT)0;
    }
    
    int rows = vsrc.rows();
    int srcstride = vsrc.row_stride(); // note that strides can be negative
    
    double sumxx = 0;
    TT sumx = (TT)0;

    TT const *psrc = vsrc.ptr();
    
    int block4 = 4*(rows/4);
    int row = 0;
    
    for(; row<block4; row += 4, psrc += 4*srcstride)
    {
        TT v = psrc[0];
        sumx += v;
        sumxx += core_real(v * core_conj(v));
        
        v = psrc[srcstride];
        sumx += v;
        sumxx += core_real(v * core_conj(v));
        
        v = psrc[2*srcstride];
        sumx += v;
        sumxx += core_real(v * core_conj(v));
        
        v = psrc[3*srcstride];
        sumx += v;
        sumxx += core_real(v * core_conj(v));
    }
    
    for(; row < rows; row++, psrc += srcstride)
    {
        TT v = psrc[0];
        sumx += v;
        sumxx += core_real(v * core_conj(v));
    }
 
    TT mean = sumx / TT(vsrc.count());
    if(pmeanrtn)
        *pmeanrtn = mean;
    
    if(vsrc.count()<2)
        return (TT)0;
    
    return TT(( sumxx - core_real(mean * core_conj(mean)) * vsrc.count() )/(vsrc.count()-1));
}

#define INST(TT) template TT im::core_block_reduce_variance(VecView<TT> const &vsrc, TT *pmeanrtn)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_sum_squares(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceSumSquares<TT>>(vsrc, (TT)0);
}

#define INST(TT) template TT im::core_block_reduce_sum_squares(VecView<TT> const &vsrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_squared_distance(VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    IM_CHECK_VALID(vsrc1);
    IM_CHECK_VALID(vsrc2);
    IM_CHECK_VECTOR_SIZES_MATCH(vsrc1, vsrc2);
    
    int rows = vsrc1.rows();
    
    int src1stride = vsrc1.row_stride(); // note that strides can be negative
    int src2stride = vsrc2.row_stride();
    
    double sumdiff = 0;

    TT const *psrc1 = vsrc1.ptr();
    TT const *psrc2 = vsrc2.ptr();
    
    int block4 = 4*(rows/4);
    int row = 0;
    
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
    }
    
    for(; row < rows; row++, psrc1 += src1stride, psrc2 += src2stride)
    {
        TT d = psrc1[0] - psrc2[0];
        sumdiff += core_real(d * core_conj(d));
    }

    
    return (TT)sumdiff;
}

#define INST(TT) template TT im::core_block_reduce_squared_distance(VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_multiply_add(VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
{
    IM_CHECK_VALID(vsrc1);
    IM_CHECK_VALID(vsrc2);
    IM_CHECK_VECTOR_SIZES_MATCH(vsrc1, vsrc2);
    
    TT sum = (TT)0;
    
    core_madd_loop(sum, vsrc1.ptr(), vsrc1.row_stride(), vsrc2.ptr(), vsrc2.row_stride(), vsrc1.rows());
    
    return sum;
}

#define INST(TT) template TT im::core_block_reduce_multiply_add(VecView<TT> const &vsrc1, VecView<TT> const &vsrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT, typename OP> TT core_priv_block_reduce_op(im::MtxView<TT> const &msrc, TT const &tally_start)
{
    IM_CHECK_VALID(msrc);
    
    // TODO: add SSE here
    
    // Fallback implementation
    TT tally = tally_start;
    
    int rows = msrc.rows();
    int cols = msrc.cols();
    
    if(rows>cols)
    {
        int srcstride = msrc.row_stride(); // note that strides can be negative
        
        for(int col = 0; col < cols; col++)
        {
            TT const *psrc = msrc.ptr(0,col);
            
            int block4 = 4*(rows/4);
            int row = 0;
            
            for(; row<block4; row += 4, psrc += 4*srcstride)
            {
                tally = OP::apply(psrc[0], tally);
                tally = OP::apply(psrc[srcstride], tally);
                tally = OP::apply(psrc[2*srcstride], tally);
                tally = OP::apply(psrc[3*srcstride], tally);
            }
            
            for(; row < rows; row++, psrc += srcstride)
                tally = OP::apply(psrc[0], tally);
        }
    }
    else
    {
        int srcstride = msrc.col_stride(); // note that strides can be negative
        
        for(int row = 0; row < rows; row++)
        {
            TT const *psrc = msrc.ptr(row,0);
            
            int block4 = 4*(cols/4);
            int col = 0;
            
            for(; col<block4; col += 4, psrc += 4*srcstride)
            {
                tally = OP::apply(psrc[0], tally);
                tally = OP::apply(psrc[srcstride], tally);
                tally = OP::apply(psrc[2*srcstride], tally);
                tally = OP::apply(psrc[3*srcstride], tally);
            }
            
            for(; col < cols; col++, psrc += srcstride)
                tally = OP::apply(psrc[0], tally);
        }
    }
    
    return tally;
}


template <typename TT> TT im::core_block_reduce_add(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceAdd<TT>>(msrc, (TT)0);
}

#define INST(TT) template TT im::core_block_reduce_add(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_multiply(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceMultiply<TT>>(msrc, (TT)1);
}

#define INST(TT) template TT im::core_block_reduce_multiply(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_mean(MtxView<TT> const &msrc)
{
    if(msrc.count()<1)
        return (TT)0;
    
    return core_priv_block_reduce_op<TT, ReduceAdd<TT>>(msrc, (TT)0) / (TT)msrc.count();
}

#define INST(TT) template TT im::core_block_reduce_mean(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_variance(MtxView<TT> const &msrcin, TT *pmeanrtn)
{
    IM_CHECK_VALID(msrcin);
    
    if(msrcin.count()<1)
    {
        if(pmeanrtn)
            *pmeanrtn = (TT)0;
        return (TT)0;
    }
    
    MtxView<TT> msrc = msrcin;
    
    if(msrc.rows()>msrc.cols())
        msrc = msrcin.t();
    
    int rows = msrc.rows();
    int cols = msrc.cols();
    int srcstride = msrc.col_stride(); // note that strides can be negative
    
    double sumxx = 0;
    TT sumx = (TT)0;
    
    for(int row = 0; row < rows; row++)
    {
        TT const *psrc = msrc.ptr(row,0);
        
        int block4 = 4*(cols/4);
        int col = 0;
        
        for(; col<block4; col += 4, psrc += 4*srcstride)
        {
            TT v = psrc[0];
            sumx += v;
            sumxx += core_real(v * core_conj(v));
            
            v = psrc[srcstride];
            sumx += v;
            sumxx += core_real(v * core_conj(v));
            
            v = psrc[2*srcstride];
            sumx += v;
            sumxx += core_real(v * core_conj(v));
            
            v = psrc[3*srcstride];
            sumx += v;
            sumxx += core_real(v * core_conj(v));
        }
        
        for(; col < cols; col++, psrc += srcstride)
        {
            TT v = psrc[0];
            sumx += v;
            sumxx += core_real(v * core_conj(v));
        }
    }
    
    TT mean = sumx / TT(msrc.count());
    if(pmeanrtn)
        *pmeanrtn = mean;
    
    if(msrcin.count()<2)
        return (TT)0;
    
    return TT(( sumxx - core_real(mean * core_conj(mean)) * msrc.count() )/(msrc.count()-1));
}

#define INST(TT) template TT im::core_block_reduce_variance(MtxView<TT> const &msrcin, TT *pmeanrtn)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_sum_squares(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceSumSquares<TT>>(msrc, (TT)0);
}

#define INST(TT) template TT im::core_block_reduce_sum_squares(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_squared_distance(MtxView<TT> const &msrc1in, MtxView<TT> const &msrc2in)
{
    IM_CHECK_VALID(msrc1in);
    IM_CHECK_VALID(msrc2in);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc1in, msrc2in);
    
    MtxView<TT> msrc1 = msrc1in;
    MtxView<TT> msrc2 = msrc2in;
    
    if(msrc1.rows()>msrc1.cols())
    {
        msrc1 = msrc1in.t();
        msrc2 = msrc2in.t();
    }
    
    int rows = msrc1.rows();
    int cols = msrc1.cols();
    
    int src1stride = msrc1.col_stride(); // note that strides can be negative
    int src2stride = msrc2.col_stride();
    
    double sumdiff = 0;
    
    for(int row = 0; row < rows; row++)
    {
        TT const *psrc1 = msrc1.ptr(row,0);
        TT const *psrc2 = msrc2.ptr(row,0);
        
        int block4 = 4*(cols/4);
        int col = 0;
        
        for(; col<block4; col += 4, psrc1 += 4*src1stride, psrc2 += 4*src2stride)
        {
            TT d = psrc1[0] - psrc2[0];
            sumdiff += core_real(d * core_conj(d));
            
            d = psrc1[src1stride] - psrc2[src2stride];
            sumdiff += core_real(d * core_conj(d));
            
            d = psrc1[2*src1stride] - psrc2[2*src2stride];
            sumdiff += core_real(d * core_conj(d));
            
            d = psrc1[3*src1stride] - psrc2[3*src2stride];
            sumdiff += core_real(d * core_conj(d));
        }
        
        for(; col < cols; col++, psrc1 += src1stride, psrc2 += src2stride)
        {
            TT d = psrc1[0] - psrc2[0];
            sumdiff += core_real(d * core_conj(d));
        }
    }
    
    return (TT)sumdiff;
}

#define INST(TT) template TT im::core_block_reduce_squared_distance(MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_multiply_add(MtxView<TT> const &msrc1in, MtxView<TT> const &msrc2in)
{
    IM_CHECK_VALID(msrc1in);
    IM_CHECK_VALID(msrc2in);
    IM_CHECK_MATRIX_SIZES_MATCH(msrc1in,msrc2in);
    
    MtxView<TT> msrc1 = msrc1in;
    MtxView<TT> msrc2 = msrc2in;
    
    if(msrc1.rows()>msrc1.cols())
    {
        msrc1 = msrc1in.t();
        msrc2 = msrc2in.t();
    }
    
    TT sum = (TT)0;
    
    for(int row = 0; row < msrc1.rows(); row++)
        core_madd_loop(sum, msrc1.ptr(row,0), msrc1.col_stride(), msrc2.ptr(row,0), msrc2.col_stride(), msrc1.cols());
    
    return sum;
}

#define INST(TT) template TT im::core_block_reduce_multiply_add(MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> TT im::core_block_reduce_max(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceMax<TT>>(vsrc, vsrc(0));
}

#define INST(TT) template TT im::core_block_reduce_max(VecView<TT> const &vsrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_min(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceMin<TT>>(vsrc, vsrc(0));
}

#define INST(TT) template TT im::core_block_reduce_min(VecView<TT> const &vsrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_max_abs(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceMaxAbs<TT>>(vsrc, std::abs(vsrc(0)));
}

#define INST(TT) template TT im::core_block_reduce_max_abs(VecView<TT> const &vsrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_min_abs(VecView<TT> const &vsrc)
{
    return core_priv_block_reduce_op<TT, ReduceMinAbs<TT>>(vsrc, std::abs(vsrc(0)));
}

#define INST(TT) template TT im::core_block_reduce_min_abs(VecView<TT> const &vsrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_median(VecView<TT> const &vsrc, int k)
{
    IM_CHECK_VALID(vsrc);
    
    int data_count = vsrc.count();
    if(data_count<1)
        return TT(0);
    
    if(k<0)
        k = data_count/2;
    if(k>=data_count)
        k = data_count - 1;
    
    std::vector<int> vindex;
    vindex.resize(data_count);
    int i;
    for(i=0; i<data_count; i++)
        vindex[i] = i;
    
    int left = 0;
    int right = data_count-1;
    
    while(true)
    {
        if(right <= left + 1)
        {
            if(right == left+1 && core_real(vsrc.index(vindex[right])) < core_real(vsrc.index(vindex[left])))
                std::swap(vindex[right], vindex[left]);
            
            return vsrc.index(vindex[k]);
        }
        else
        {
            std::swap(vindex[(left + right)/2], vindex[left+1]);
            
            if(core_real(vsrc.index(vindex[left])) > core_real(vsrc.index(vindex[right])))
                std::swap(vindex[left], vindex[right]);
            
            if(core_real(vsrc.index(vindex[left+1])) > core_real(vsrc.index(vindex[right])))
                std::swap(vindex[left+1], vindex[right]);
            
            if(core_real(vsrc.index(vindex[left])) > core_real(vsrc.index(vindex[left+1])))
                std::swap(vindex[left], vindex[left+1]);
            
            int i = left+1;
            int j = right;
            TT z = vsrc.index(vindex[left+1]);
            
            while (true)
            {
                do i++; while(core_real(vsrc.index(vindex[i]))<core_real(z));
                do j--; while(core_real(vsrc.index(vindex[j]))>core_real(z));
                if(j<i)
                    break;
                std::swap(vindex[i], vindex[j]);
            }
            
            std::swap(vindex[left+1], vindex[j]);
            
            if(j>=k)
                right = j-1;
            if(j<=k)
                left = i;
        }
    }
    
    return (TT)0;
}

#define INST(TT) template TT im::core_block_reduce_median(VecView<TT> const &vsrc, int k)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST


template <typename TT> TT im::core_block_reduce_max(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceMax<TT>>(msrc, msrc(0,0));
}

#define INST(TT) template TT im::core_block_reduce_max(MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_min(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceMin<TT>>(msrc, msrc(0,0));
}

#define INST(TT) template TT im::core_block_reduce_min(MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_max_abs(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceMaxAbs<TT>>(msrc, std::abs(msrc(0,0)));
}

#define INST(TT) template TT im::core_block_reduce_max_abs(MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_min_abs(MtxView<TT> const &msrc)
{
    return core_priv_block_reduce_op<TT, ReduceMinAbs<TT>>(msrc, std::abs(msrc(0,0)));
}

#define INST(TT) template TT im::core_block_reduce_min_abs(MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> TT im::core_block_reduce_median(MtxView<TT> const &msrc, int k)
{
    IM_CHECK_VALID(msrc);
    
    int data_count = msrc.count();
    if(data_count<1)
        return TT(0);
    
    if(k<0)
        k = data_count/2;
    if(k>=data_count)
        k = data_count - 1;
    
    std::vector<int> vindex;
    vindex.resize(data_count);
    int i;
    for(i=0; i<data_count; i++)
        vindex[i] = i;
    
    int left = 0;
    int right = data_count-1;
    
    while(true)
    {
        if(right <= left + 1)
        {
            if(right == left+1 && core_real(msrc.index(vindex[right])) < core_real(msrc.index(vindex[left])))
                std::swap(vindex[right], vindex[left]);
            
            return msrc.index(vindex[k]);
        }
        else
        {
            std::swap(vindex[(left + right)/2], vindex[left+1]);
            
            if(core_real(msrc.index(vindex[left])) > core_real(msrc.index(vindex[right])))
                std::swap(vindex[left], vindex[right]);
            
            if(core_real(msrc.index(vindex[left+1])) > core_real(msrc.index(vindex[right])))
                std::swap(vindex[left+1], vindex[right]);
            
            if(core_real(msrc.index(vindex[left])) > core_real(msrc.index(vindex[left+1])))
                std::swap(vindex[left], vindex[left+1]);
            
            int i = left+1;
            int j = right;
            TT z = msrc.index(vindex[left+1]);
            
            while (true)
            {
                do i++; while(core_real(msrc.index(vindex[i]))<core_real(z));
                do j--; while(core_real(msrc.index(vindex[j]))>core_real(z));
                if(j<i)
                    break;
                std::swap(vindex[i], vindex[j]);
            }
            
            std::swap(vindex[left+1], vindex[j]);
            
            if(j>=k)
                right = j-1;
            if(j<=k)
                left = i;
        }
    }
    
    return (TT)0;
}

#define INST(TT) template TT im::core_block_reduce_median(MtxView<TT> const &msrc, int k)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST


template <typename TT> void im::core_block_reduce_rows_add(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_add(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_add(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_multiply(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_multiply(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_multiply(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_max(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_max(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_max(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_min(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_min(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_min(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_max_abs(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_max_abs(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_max_abs(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_min_abs(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_min_abs(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_min_abs(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_median(VecView<TT> vdst, MtxView<TT> const &msrc, int k)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_median(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_median(VecView<TT> vdst, MtxView<TT> const &msrc, int k)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_mean(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_mean(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_mean(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_variance(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_variance(msrc.row(row), (TT *)NULL);
}

#define INST(TT) template void im::core_block_reduce_rows_variance(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_sum_squares(VecView<TT> vdst, MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_sum_squares(msrc.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_sum_squares(VecView<TT> vdst, MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_squared_distance(VecView<TT> vdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc1.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_squared_distance(msrc1.row(row), msrc2.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_squared_distance(VecView<TT> vdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_reduce_rows_multiply_add(VecView<TT> vdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
{
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vdst.rows()==msrc1.rows());
    for(int row = 0; row<vdst.rows(); row++)
        vdst(row) = core_block_reduce_multiply_add(msrc1.row(row), msrc2.row(row));
}

#define INST(TT) template void im::core_block_reduce_rows_multiply_add(VecView<TT> vdst, MtxView<TT> const &msrc1, MtxView<TT> const &msrc2)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

