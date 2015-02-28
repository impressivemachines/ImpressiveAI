//
//  block_blas1.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/9/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT> TT im::core_block_blas_asum(VecView<TT> const &vx)
{
    return core_block_reduce_add(vx);
}

#define INST(TT) template TT im::core_block_blas_asum(VecView<TT> const &vx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_axpy(VecView<TT> vy, VecView<TT> const &vx, const TT &alpha)
{
    IM_CHECK_VALID(vx);
    IM_CHECK_VALID(vy);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    if(alpha==(TT)0)
        return;
    
    core_sadd_loop(vy.ptr(), vy.row_stride(), vx.ptr(), vx.row_stride(), alpha, vx.rows());
}

#define INST(TT) template void im::core_block_blas_axpy(VecView<TT> vy, VecView<TT> const &vx, const TT &alpha)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_copy(VecView<TT> vy, VecView<TT> const &vx)
{
    core_block_copy(vy,vx);
}

#define INST(TT) template void im::core_block_blas_copy(VecView<TT> vy, VecView<TT> const &vx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> TT im::core_block_blas_dot(VecView<TT> const &vx, VecView<TT> const &vy)
{
    return core_block_reduce_multiply_add(vx,vy);
}

#define INST(TT) template TT im::core_block_blas_dot(VecView<TT> const &vx, VecView<TT> const &vy)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> TT im::core_block_blas_dotc(VecView<TT> const &vx, VecView<TT> const &vy)
{
    IM_CHECK_VALID(vx);
    IM_CHECK_VALID(vy);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    TT sum = (TT)0;
    
    core_maddconj_loop(sum, vx.ptr(), vx.row_stride(), vy.ptr(), vy.row_stride(), vx.rows());
    
    return sum;
}

#define INST(TT) template TT im::core_block_blas_dotc(VecView<TT> const &vx, VecView<TT> const &vy)
INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

double im::core_block_blas_sdot(VecView<float> const &vx, VecView<float> const &vy)
{
    IM_CHECK_VALID(vx);
    IM_CHECK_VALID(vy);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    double sum = 0.0;
    
    core_madd_loop(sum, vx.ptr(), vx.row_stride(), vy.ptr(), vy.row_stride(), vx.rows());
    
    return sum;
}

//

template <typename TT> TT im::core_block_blas_nrm2(VecView<TT> const &vx)
{
    IM_CHECK_VALID(vx);
    
    TT scale = (TT)0;
    TT ssq = (TT)1;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT x = vx(i);
        if(x != (TT)0)
        {
            TT ax = std::abs(x);
            if(scale<ax)
            {
                TT sax = scale / ax;
                ssq = (TT)1 + ssq * sax * sax;
                scale = ax;
            }
            else
            {
                TT axs = ax / scale;
                ssq += axs * axs;
            }
        }
    }
    
    return scale * std::sqrt(ssq);
}

template <typename TT> std::complex<TT> im::core_block_blas_nrm2(VecView<std::complex<TT>> const &vx)
{
    IM_CHECK_VALID(vx);
    
    TT scale = (TT)0;
    TT ssq = (TT)1;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT x = vx(i).real();
        if(x != (TT)0)
        {
            TT ax = std::abs(x);
            if(scale<ax)
            {
                TT sax = scale / ax;
                ssq = (TT)1 + ssq * sax * sax;
                scale = ax;
            }
            else
            {
                TT axs = ax / scale;
                ssq += axs * axs;
            }
        }
        
        TT y = vx(i).imag();
        if(y != (TT)0)
        {
            TT ay = std::abs(y);
            if(scale<ay)
            {
                TT say = scale / ay;
                ssq = (TT)1 + ssq * say * say;
                scale = ay;
            }
            else
            {
                TT ays = ay / scale;
                ssq += ays * ays;
            }
        }
    }
    
    return std::complex<TT>(scale * std::sqrt(ssq));
}

#define INST(TT) template TT im::core_block_blas_nrm2(VecView<TT> const &vx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> TT im::core_block_blas_nrm2(MtxView<TT> const &mx)
{
    IM_CHECK_VALID(mx);
    
    TT scale = (TT)0;
    TT ssq = (TT)1;
    
    for(int i=0; i<mx.rows(); i++)
        for(int j=0; j<mx.cols(); j++)
        {
            TT x = mx(i,j);
            if(x != (TT)0)
            {
                TT ax = std::abs(x);
                if(scale<ax)
                {
                    TT sax = scale / ax;
                    ssq = (TT)1 + ssq * sax * sax;
                    scale = ax;
                }
                else
                {
                    TT axs = ax / scale;
                    ssq += axs * axs;
                }
            }
        }
    
    return scale * std::sqrt(ssq);
}

template <typename TT> std::complex<TT> im::core_block_blas_nrm2(MtxView<std::complex<TT>> const &mx)
{
    IM_CHECK_VALID(mx);
    
    TT scale = (TT)0;
    TT ssq = (TT)1;
    
    for(int i=0; i<mx.rows(); i++)
        for(int j=0; j<mx.cols(); j++)
        {
            TT x = mx(i,j).real();
            if(x != (TT)0)
            {
                TT ax = std::abs(x);
                if(scale<ax)
                {
                    TT sax = scale / ax;
                    ssq = (TT)1 + ssq * sax * sax;
                    scale = ax;
                }
                else
                {
                    TT axs = ax / scale;
                    ssq += axs * axs;
                }
            }
            
            TT y = mx(i,j).imag();
            if(y != (TT)0)
            {
                TT ay = std::abs(y);
                if(scale<ay)
                {
                    TT say = scale / ay;
                    ssq = (TT)1 + ssq * say * say;
                    scale = ay;
                }
                else
                {
                    TT ays = ay / scale;
                    ssq += ays * ays;
                }
            }
        }
    
    return std::complex<TT>(scale * std::sqrt(ssq));
}

#define INST(TT) template TT im::core_block_blas_nrm2(MtxView<TT> const &mx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_rot(VecView<TT> vx, VecView<TT> vy, TT const &c, TT const &s)
{
    IM_CHECK_VALID(vx);
    IM_CHECK_VALID(vy);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT x = vx(i);
        TT y = vy(i);
        vx(i) = c * x + s * y;
        vy(i) = -s * x + c * y;
    }
}

#define INST(TT) template void im::core_block_blas_rot(VecView<TT> vx, VecView<TT> vy, TT const &c, TT const &s)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_rotg(TT &a, TT &b, TT &c, TT &s)
{
    TT roe = (std::abs(a) > std::abs(b) ? a : b);
    TT scale = std::abs(a) + std::abs(b);
    
    TT r, z;
    
    if(scale != (TT)0)
    {
        TT aos = a / scale;
        TT bos = b / scale;
        r = scale * std::sqrt(aos * aos + bos * bos);
        r = core_sgn(roe) * r;
        c = a / r;
        s = b / r;
        z = (TT)1;
        if(std::abs(a) > std::abs(b))
            z = s;
        if(std::abs(b) >= std::abs(a) && c != (TT)0)
            z = (TT)1 / c;
    }
    else
    {
        c = (TT)1;
        s = (TT)0;
        r = (TT)0;
        z = (TT)0;
    }
    
    a = r;
    b = z;
}

#define INST(TT) template void im::core_block_blas_rotg(TT &a, TT &b, TT &c, TT &s)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_rotm(VecView<TT> vx, VecView<TT> vy, TT const *pparam)
{
    IM_CHECK_VALID(vx);
    IM_CHECK_VALID(vy);
    IM_CHECK_VECTOR_SIZES_MATCH(vx, vy);
    
    TT h11, h21, h12, h22;
    
    if(pparam[0] == (TT)(-1))
    {
        h11 = pparam[1];
        h21 = pparam[2];
        h12 = pparam[3];
        h22 = pparam[4];
    }
    else if(pparam[0] == (TT)(0))
    {
        h11 = (TT)1;
        h21 = pparam[2];
        h12 = pparam[3];
        h22 = (TT)1;
    }
    else if(pparam[0] == (TT)(1))
    {
        h11 = pparam[1];
        h21 = (TT)(-1);
        h12 = (TT)(1);
        h22 = pparam[4];
    }
    else if(pparam[0] == (TT)(2))
    {
        return;
    }
    else
        IM_THROW_ARGUMENTS;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT w = vx(i);
        TT z = vy(i);
        vx(i) = h11 * w + h12 * z;
        vy(i) = h21 * w + h22 * z;
    }
}

#define INST(TT) template void im::core_block_blas_rotm(VecView<TT> vx, VecView<TT> vy, TT const *pparam)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_rotmg(TT *pparam, TT &d1, TT &d2, TT &b1, TT const &b2)
{
    TT const G = 4096;
    TT G2 = G * G;
    TT D1 = d1;
    TT D2 = d2;
    TT x = b1;
    TT y = b2;
    TT h11, h12, h21, h22, u;
    TT c, s;

    if(D1 < 0)
    {
        pparam[0] = -1;
        pparam[1] = 0;
        pparam[2] = 0;
        pparam[3] = 0;
        pparam[4] = 0;
        d1 = 0;
        d2 = 0;
        b1 = 0;
        return;
    }
    
    if(D2 * y == 0)
    {
        pparam[0] = -2;
        return;
    }
    
    c = std::abs(D1 * x * x);
    s = std::abs(D2 * y * y);
    
    if(c > s)
    {
        pparam[0] = 0;
        
        h11 = 1;
        h12 = (D2 * y) / (D1 * x);
        h21 = -y / x;
        h22 = 1;
        
        u = 1 - h21 * h12;
        
        if(u <= 0)
        {
            pparam[0] = -1;
            pparam[1] = 0;
            pparam[2] = 0;
            pparam[3] = 0;
            pparam[4] = 0;
            d1 = 0;
            d2 = 0;
            b1 = 0;
            return;
        }
        
        D1 /= u;
        D2 /= u;
        x *= u;
    }
    else
    {
        if(D2 * y * y < 0)
        {
            pparam[0] = -1;
            pparam[1] = 0;
            pparam[2] = 0;
            pparam[3] = 0;
            pparam[4] = 0;
            d1 = 0;
            d2 = 0;
            b1 = 0;
            return;
        }
        
        pparam[0] = 1;
        
        h11 = (D1 * x) / (D2 * y);
        h12 = 1;
        h21 = -1;
        h22 = x / y;
        
        u = 1 + h11 * h22;
        
        D1 /= u;
        D2 /= u;
        
        {
            TT tmp = D2;
            D2 = D1;
            D1 = tmp;
        }
        
        x = y * u;
    }

    while (D1 <= (TT)1.0 / G2 && D1 != 0)
    {
        pparam[0] = -1;
        D1 *= G2;
        x /= G;
        h11 /= G;
        h12 /= G;
    }
    
    while (D1 >= G2)
    {
        pparam[0] = -1;
        D1 /= G2;
        x *= G;
        h11 *= G;
        h12 *= G;
    }

    while (std::abs(D2) <= (TT)1.0 / G2 && D2 != 0)
    {
        pparam[0] = -1;
        D2 *= G2;
        h21 /= G;
        h22 /= G;
    }
    
    while (std::abs(D2) >= G2)
    {
        pparam[0] = -1;
        D2 /= G2;
        h21 *= G;
        h22 *= G;
    }
    
    d1 = D1;
    d2 = D2;
    b1 = x;
    
    if(pparam[0] == -1)
    {
        pparam[1] = h11;
        pparam[2] = h21;
        pparam[3] = h12;
        pparam[4] = h22;
    }
    else if(pparam[0] == 0)
    {
        pparam[2] = h21;
        pparam[3] = h12;
    }
    else if(pparam[0] == 1)
    {
        pparam[1] = h11;
        pparam[4] = h22;
    }
}

#define INST(TT) template void im::core_block_blas_rotmg(TT *pparam, TT &d1, TT &d2, TT &x, TT const &y)
INST(float); INST(double);
#undef INST

//

template <typename TT> void im::core_block_blas_scal(VecView<TT> vx, TT const &alpha)
{
    core_block_scale(vx,vx,alpha);
}

#define INST(TT) template void im::core_block_blas_scal(VecView<TT> vx, TT const &alpha)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void im::core_block_blas_swap(VecView<TT> vx, VecView<TT> vy)
{
    core_block_exchange(vx,vy);
}

#define INST(TT) template void im::core_block_blas_swap(VecView<TT> vx, VecView<TT> vy)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> int im::core_block_blas_iamax(VecView<TT> const &vx)
{
    IM_CHECK_VALID(vx);
    
    TT maxval = (TT)0;
    int maxind = 0;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT const a = std::abs(vx(i));
        if(a > maxval)
        {
            maxval = a;
            maxind = i;
        }
    }
    
    return maxind;
}

template <typename TT> int im::core_block_blas_iamax(VecView<std::complex<TT>> const &vx)
{
    IM_CHECK_VALID(vx);
    
    TT maxval = (TT)0;
    int maxind = 0;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT const a = core_sumabs(vx(i));
        if(a > maxval)
        {
            maxval = a;
            maxind = i;
        }
    }
    
    return maxind;
}

#define INST(TT) template int im::core_block_blas_iamax(VecView<TT> const &vx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> int im::core_block_blas_iamin(VecView<TT> const &vx)
{
    IM_CHECK_VALID(vx);
    
    TT minval = (TT)0;
    int minind = 0;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT const a = std::abs(vx(i));
        if(a < minval)
        {
            minval = a;
            minind = i;
        }
    }
    
    return minind;
}

template <typename TT> int im::core_block_blas_iamin(VecView<std::complex<TT>> const &vx)
{
    IM_CHECK_VALID(vx);
    
    TT minval = (TT)0;
    int minind = 0;
    
    for(int i=0; i<vx.rows(); i++)
    {
        TT const a = core_sumabs(vx(i));
        if(a < minval)
        {
            minval = a;
            minind = i;
        }
    }
    
    return minind;
}

#define INST(TT) template int im::core_block_blas_iamin(VecView<TT> const &vx)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST



