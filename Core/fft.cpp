//
//  fft.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/26/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

#define SPECIALCASE_POWEROF2    1

template <typename TT> void im::FFTComplex<TT>::init(int n)
{
    m_factors.clear();
    m_offsets.clear();
    m_n = 0;
    
    IM_CHECK_ARGS(n>=0);

    if(n<2)
        return;
    
    m_n = n;
    
#if SPECIALCASE_POWEROF2
    if(core_is_pow2(m_n))
        return;
#endif
    
    m_twiddle.resize(m_n);
    m_buffer.resize(m_n);
    
    // special case factors
    factor_reduce(n, 5);
    factor_reduce(n, 4);
    factor_reduce(n, 3);
    factor_reduce(n, 2);
    
    // general odd factors
    for(int factor = 7; n>1; factor+=2)
        factor_reduce(n, factor);
    
    m_offsets.resize(m_factors.size());
    
    int offset = 0;
    int product = 1;
    double dtheta = -2.0*CONST_PI / m_n;
    
    // generate all the trig twiddle factors
    for(int i=0; i<m_factors.size(); i++)
    {
        int factor = m_factors[i];
        m_offsets[i] = offset;
        int product_curr = product;
        product *= factor;
        int trigcount = m_n / product;
        
        for(int j=1; j<factor; j++)
        {
            int a = 0;
            for(int k=1; k<=trigcount; k++)
            {
                a += j * product_curr;
                a = a % m_n;
                
                double theta = dtheta * a;
                m_twiddle[offset].real((TT)std::cos(theta));
                m_twiddle[offset].imag((TT)std::sin(theta));
                offset++;
            }
        }
    }
    
    IM_CHECK(offset<=m_n);
}

template <typename TT> void im::FFTComplex<TT>::factor_reduce(int &n, int factor)
{
    while(factor!=1 && (n % factor)==0)
    {
        m_factors.push_back(factor);
        n /= factor;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform_radix_2(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase)
{
    int const factor = 2;
    int const trigcount = m_n / product;
    int const offset_j = product / factor;
    int const offset_i = m_n / factor;
    int const jinc = (factor-1) * offset_j;
    std::complex<TT> const *ptwid1 = m_twiddle.data() + twiddlebase;
    int const sign = (int)m_direction;
    int const msign = -sign;
    
    int i = 0; // input index
    int j = 0; // output index
    
    // FFT2 is matrix multiply by [1 1; 1 -1]
    // followed by twiddle factor scaling [1 W1]
    
    for(int k=0; k<trigcount; k++)
    {
        // W twiddle factor
        TT w_re, w_im;
        
        if(k==0)
        {
            w_re = 1;
            w_im = 0;
        }
        else
        {
            w_re = ptwid1[k-1].real();
            w_im = msign * ptwid1[k-1].imag();
        }
        
        for(int u=0; u<offset_j; u++, i++, j++)
        {
            TT const z0_re = vsrc(i).real();
            TT const z0_im = vsrc(i).imag();
            TT const z1_re = vsrc(i+offset_i).real();
            TT const z1_im = vsrc(i+offset_i).imag();
            
            TT const x0_re = z0_re + z1_re;
            TT const x0_im = z0_im + z1_im;
            TT const x1_re = z0_re - z1_re;
            TT const x1_im = z0_im - z1_im;
            
            // 1 * X0
            vdst(j).real(x0_re);
            vdst(j).imag(x0_im);
            
            // W * X1
            vdst(j+offset_j).real(w_re * x1_re - w_im * x1_im);
            vdst(j+offset_j).imag(w_re * x1_im + w_im * x1_re);
        }
        
        j += jinc;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform_radix_3(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase)
{
    int const factor = 3;
    int const trigcount = m_n / product;
    int const offset_i = m_n / factor;
    int const offset_j = product / factor;
    int const jinc = (factor-1) * offset_j;
    std::complex<TT> const *ptwid1 = m_twiddle.data() + twiddlebase;
    std::complex<TT> const *ptwid2 = m_twiddle.data() + twiddlebase + trigcount;
    int const sign = (int)m_direction;
    int const msign = -sign;
    
    TT const sqrt3_2 = (TT)CONST_SQRT3_2;
    
    int i = 0; // input index
    int j = 0; // output index
    
    // FFT3 is matrix multiply by 3x3 coef matrix
    // followed by twiddle factor scaling [1 W1 W2]
    
    for(int k=0; k<trigcount; k++)
    {
        // W twiddle factor
        TT w1_re, w1_im, w2_re, w2_im;
        
        if(k==0)
        {
            w1_re = 1;
            w1_im = 0;
            w2_re = 1;
            w2_im = 0;
        }
        else
        {
            w1_re = ptwid1[k-1].real();
            w1_im = msign * ptwid1[k-1].imag();
            w2_re = ptwid2[k-1].real();
            w2_im = msign * ptwid2[k-1].imag();
        }
        
        for(int u=0; u<offset_j; u++, i++, j++)
        {
            TT const z0_re = vsrc(i).real();
            TT const z0_im = vsrc(i).imag();
            TT const z1_re = vsrc(i+offset_i).real();
            TT const z1_im = vsrc(i+offset_i).imag();
            TT const z2_re = vsrc(i+2*offset_i).real();
            TT const z2_im = vsrc(i+2*offset_i).imag();
            
            // X = W3x3 * Z
            TT const t1_re = z1_re + z2_re;
            TT const t1_im = z1_im + z2_im;
            TT const t2_re = z0_re - t1_re / 2.0;
            TT const t2_im = z0_im - t1_im / 2.0;
            TT const t3_re = sign * sqrt3_2 * (z1_re - z2_re);
            TT const t3_im = sign * sqrt3_2 * (z1_im - z2_im);
            
            TT const x0_re = z0_re + t1_re;
            TT const x0_im = z0_im + t1_im;
            TT const x1_re = t2_re - t3_im;
            TT const x1_im = t2_im + t3_re;
            TT const x2_re = t2_re + t3_im;
            TT const x2_im = t2_im - t3_re;
            
            // 1 * X0
            vdst(j).real(x0_re);
            vdst(j).imag(x0_im);
            
            // W1 * X1
            vdst(j+offset_j).real(w1_re * x1_re - w1_im * x1_im);
            vdst(j+offset_j).imag(w1_re * x1_im + w1_im * x1_re);
            
            // W2 * X2
            vdst(j+2*offset_j).real(w2_re * x2_re - w2_im * x2_im);
            vdst(j+2*offset_j).imag(w2_re * x2_im + w2_im * x2_re);
        }
        
        j += jinc;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform_radix_4(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase)
{
    int const factor = 4;
    int const trigcount = m_n / product;
    int const offset_i = m_n / factor;
    int const offset_j = product / factor;
    int const jinc = (factor-1) * offset_j;
    std::complex<TT> const *ptwid1 = m_twiddle.data() + twiddlebase;
    std::complex<TT> const *ptwid2 = m_twiddle.data() + twiddlebase + trigcount;
    std::complex<TT> const *ptwid3 = m_twiddle.data() + twiddlebase + 2*trigcount;
    int const sign = (int)m_direction;
    int const msign = -sign;
    
    int i = 0; // input index
    int j = 0; // output index
    
    // FFT4 is matrix multiply by 4x4 coef matrix
    // followed by twiddle factor scaling [1 W1 W2 W3]
    
    for(int k=0; k<trigcount; k++)
    {
        // W twiddle factor
        TT w1_re, w1_im, w2_re, w2_im, w3_re, w3_im;
        
        if(k==0)
        {
            w1_re = 1;
            w1_im = 0;
            w2_re = 1;
            w2_im = 0;
            w3_re = 1;
            w3_im = 0;
        }
        else
        {
            w1_re = ptwid1[k-1].real();
            w1_im = msign * ptwid1[k-1].imag();
            w2_re = ptwid2[k-1].real();
            w2_im = msign * ptwid2[k-1].imag();
            w3_re = ptwid3[k-1].real();
            w3_im = msign * ptwid3[k-1].imag();
        }
        
        for(int u=0; u<offset_j; u++, i++, j++)
        {
            TT const z0_re = vsrc(i).real();
            TT const z0_im = vsrc(i).imag();
            TT const z1_re = vsrc(i+offset_i).real();
            TT const z1_im = vsrc(i+offset_i).imag();
            TT const z2_re = vsrc(i+2*offset_i).real();
            TT const z2_im = vsrc(i+2*offset_i).imag();
            TT const z3_re = vsrc(i+3*offset_i).real();
            TT const z3_im = vsrc(i+3*offset_i).imag();
            
            // X = W4x4 * Z
            TT const t1_re = z0_re + z2_re;
            TT const t1_im = z0_im + z2_im;
            TT const t2_re = z1_re + z3_re;
            TT const t2_im = z1_im + z3_im;
            TT const t3_re = z0_re - z2_re;
            TT const t3_im = z0_im - z2_im;
            TT const t4_re = sign * (z1_re - z3_re);
            TT const t4_im = sign * (z1_im - z3_im);
            
            TT const x0_re = t1_re + t2_re;
            TT const x0_im = t1_im + t2_im;
            TT const x1_re = t3_re - t4_im;
            TT const x1_im = t3_im + t4_re;
            TT const x2_re = t1_re - t2_re;
            TT const x2_im = t1_im - t2_im;
            TT const x3_re = t3_re + t4_im;
            TT const x3_im = t3_im - t4_re;
            
            // 1 * X0
            vdst(j).real(x0_re);
            vdst(j).imag(x0_im);
            
            // W1 * X1
            vdst(j+offset_j).real(w1_re * x1_re - w1_im * x1_im);
            vdst(j+offset_j).imag(w1_re * x1_im + w1_im * x1_re);
            
            // W2 * X2
            vdst(j+2*offset_j).real(w2_re * x2_re - w2_im * x2_im);
            vdst(j+2*offset_j).imag(w2_re * x2_im + w2_im * x2_re);
            
            // W3 * X3
            vdst(j+3*offset_j).real(w3_re * x3_re - w3_im * x3_im);
            vdst(j+3*offset_j).imag(w3_re * x3_im + w3_im * x3_re);
        }
        
        j += jinc;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform_radix_5(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase)
{
    int const factor = 5;
    int const trigcount = m_n / product;
    int const offset_i = m_n / factor;
    int const offset_j = product / factor;
    int const jinc = (factor-1) * offset_j;
    std::complex<TT> const *ptwid1 = m_twiddle.data() + twiddlebase;
    std::complex<TT> const *ptwid2 = m_twiddle.data() + twiddlebase + trigcount;
    std::complex<TT> const *ptwid3 = m_twiddle.data() + twiddlebase + 2*trigcount;
    std::complex<TT> const *ptwid4 = m_twiddle.data() + twiddlebase + 3*trigcount;
    int const sign = (int)m_direction;
    int const msign = -sign;
    
    int i = 0; // input index
    int j = 0; // output index
    
    TT const sin_2pi_5 = (TT)std::sin(2.0 * M_PI / 5.0);
    TT const sin_2pi_10 = (TT)std::sin(2.0 * M_PI / 10.0);
    TT const sqrt_5_4 = (TT)(std::sqrt(5.0)/4.0);
    
    // FFT5 is matrix multiply by 5x5 coef matrix
    // followed by twiddle factor scaling [1 W1 W2 W3 W4]
    
    for(int k=0; k<trigcount; k++)
    {
        // W twiddle factor
        TT w1_re, w1_im, w2_re, w2_im, w3_re, w3_im, w4_re, w4_im;
        
        if(k==0)
        {
            w1_re = 1;
            w1_im = 0;
            w2_re = 1;
            w2_im = 0;
            w3_re = 1;
            w3_im = 0;
            w4_re = 1;
            w4_im = 0;
        }
        else
        {
            w1_re = ptwid1[k-1].real();
            w1_im = msign * ptwid1[k-1].imag();
            w2_re = ptwid2[k-1].real();
            w2_im = msign * ptwid2[k-1].imag();
            w3_re = ptwid3[k-1].real();
            w3_im = msign * ptwid3[k-1].imag();
            w4_re = ptwid4[k-1].real();
            w4_im = msign * ptwid4[k-1].imag();
        }
        
        for(int u=0; u<offset_j; u++, i++, j++)
        {
            TT const z0_re = vsrc(i).real();
            TT const z0_im = vsrc(i).imag();
            TT const z1_re = vsrc(i+offset_i).real();
            TT const z1_im = vsrc(i+offset_i).imag();
            TT const z2_re = vsrc(i+2*offset_i).real();
            TT const z2_im = vsrc(i+2*offset_i).imag();
            TT const z3_re = vsrc(i+3*offset_i).real();
            TT const z3_im = vsrc(i+3*offset_i).imag();
            TT const z4_re = vsrc(i+4*offset_i).real();
            TT const z4_im = vsrc(i+4*offset_i).imag();
            
            // X = W5x5 * Z
            TT const t1_re = z1_re + z4_re;
            TT const t1_im = z1_im + z4_im;
            TT const t2_re = z2_re + z3_re;
            TT const t2_im = z2_im + z3_im;
            TT const t3_re = z1_re - z4_re;
            TT const t3_im = z1_im - z4_im;
            TT const t4_re = z2_re - z3_re;
            TT const t4_im = z2_im - z3_im;
            TT const t5_re = t1_re + t2_re;
            TT const t5_im = t1_im + t2_im;
            TT const t6_re = sqrt_5_4 * (t1_re - t2_re);
            TT const t6_im = sqrt_5_4 * (t1_im - t2_im);
            TT const t7_re = z0_re - t5_re / 4.0;
            TT const t7_im = z0_im - t5_im / 4.0;
            TT const t8_re = t7_re + t6_re;
            TT const t8_im = t7_im + t6_im;
            TT const t9_re = t7_re - t6_re;
            TT const t9_im = t7_im - t6_im;
            TT const t10_re = sign * (sin_2pi_5 * t3_re + sin_2pi_10 * t4_re);
            TT const t10_im = sign * (sin_2pi_5 * t3_im + sin_2pi_10 * t4_im);
            TT const t11_re = sign * (sin_2pi_10 * t3_re - sin_2pi_5 * t4_re);
            TT const t11_im = sign * (sin_2pi_10 * t3_im - sin_2pi_5 * t4_im);
            
            TT const x0_re = z0_re + t5_re;
            TT const x0_im = z0_im + t5_im;
            TT const x1_re = t8_re - t10_im;
            TT const x1_im = t8_im + t10_re;
            TT const x2_re = t9_re - t11_im;
            TT const x2_im = t9_im + t11_re;
            TT const x3_re = t9_re + t11_im;
            TT const x3_im = t9_im - t11_re;
            TT const x4_re = t8_re + t10_im;
            TT const x4_im = t8_im - t10_re;
            
            // 1 * X0
            vdst(j).real(x0_re);
            vdst(j).imag(x0_im);
            
            // W1 * X1
            vdst(j+offset_j).real(w1_re * x1_re - w1_im * x1_im);
            vdst(j+offset_j).imag(w1_re * x1_im + w1_im * x1_re);
            
            // W2 * X2
            vdst(j+2*offset_j).real(w2_re * x2_re - w2_im * x2_im);
            vdst(j+2*offset_j).imag(w2_re * x2_im + w2_im * x2_re);
            
            // W3 * X3
            vdst(j+3*offset_j).real(w3_re * x3_re - w3_im * x3_im);
            vdst(j+3*offset_j).imag(w3_re * x3_im + w3_im * x3_re);
            
            // W4 * X4
            vdst(j+4*offset_j).real(w4_re * x4_re - w4_im * x4_im);
            vdst(j+4*offset_j).imag(w4_re * x4_im + w4_im * x4_re);
        }
        
        j += jinc;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform_radix_N(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase, int factor)
{
    int const trigcount = m_n / product;
    int const offset_i = m_n / factor;
    int const offset_j = product / factor;
    int const jinc = (factor-1) * offset_j;
    int const halffactor = (factor - 1)/2 + 1;
    std::complex<TT> const  *ptwid1 = m_twiddle.data() + twiddlebase;
    
    int const sign = (int)m_direction;
    int const msign = -sign;
    
    int i; // input index
    int j; // output index
    
    int e, e1, k, k1;
    
    for(i = 0; i < offset_i; i++)
    {
        vdst(i) = vsrc(i);
    }
    
    for(e = 1; e < halffactor; e++)
    {
        for(i = 0; i < offset_i; i++)
        {
            int const u = i + e * offset_i;
            int const v = i + (factor - e) * offset_i;
            vdst(u) = vsrc(u) + vsrc(v);
            vdst(v) = vsrc(u) - vsrc(v);
        }
    }
    
    for(i=0; i<offset_i; i++)
    {
        vsrc(i) = vdst(i);
    }
    
    for(e1 = 1; e1 < halffactor; e1++)
    {
        for(i = 0; i < offset_i; i++)
        {
            int const u = i + e1 * offset_i;
            vsrc(i) += vdst(u);
        }
    }
    
    for(e = 1; e < halffactor; e++)
    {
        int idx = e * trigcount;
        int const idx_step = e * trigcount;
        int const em = e * offset_i;
        int const ecm = (factor - e) * offset_i;
        
        for(i = 0; i < offset_i; i++)
        {
            vsrc(i+em) = vdst(i);
            vsrc(i+ecm) = vdst(i);
        }
        
        for(e1 = 1; e1 < halffactor; e1++)
        {
            TT w_re, w_im;
            
            if(idx == 0) {
                w_re = 1;
                w_im = 0;
            } else {
                w_re = ptwid1[idx-1].real();
                w_im = msign * ptwid1[idx-1].imag();
            }
            
            for(i = 0; i < offset_i; i++)
            {
                int const u = i + e1 * offset_i;
                int const v = i + (factor - e1) * offset_i;
                
                TT const xp_re = vdst(u).real();
                TT const xp_im = vdst(u).imag();
                TT const xm_re = vdst(v).real();
                TT const xm_im = vdst(v).imag();
                
                TT const ap = w_re * xp_re;
                TT const am = w_im * xm_im;
                TT const bp = w_re * xp_im;
                TT const bm = w_im * xm_re;
                
                std::complex<TT> sum(ap - am, bp + bm);
                std::complex<TT> sumc(ap + am, bp - bm);
                
                vsrc(i+em) += sum;
                vsrc(i+ecm) += sumc;
            }
            idx += idx_step;
            idx %= factor * trigcount;
        }
    }
    
    i = 0;
    j = 0;
    
    for(k1 = 0; k1 < offset_j; k1++)
    {
        vdst(k1) = vsrc(k1);
    }
    
    for(e1 = 1; e1 < factor; e1++)
    {
        for(k1 = 0; k1 < offset_j; k1++)
        {
            int const u = k1 + e1 * offset_j;
            int const v = k1 + e1 * offset_i;
            vdst(u) = vsrc(v);
        }
    }
    
    i = offset_j;
    j = product;
    
    for(k = 1; k < trigcount; k++)
    {
        for(k1 = 0; k1 < offset_j; k1++)
        {
            vdst(j) = vsrc(i);
            i++;
            j++;
        }
        j += jinc;
    }
    
    i = offset_j;
    j = product;
    
    for(k = 1; k < trigcount; k++)
    {
        for(k1 = 0; k1 < offset_j; k1++)
        {
            for(e1 = 1; e1 < factor; e1++)
            {
                int const u = (e1-1)*trigcount + k-1;
                int const vi = i + e1 * offset_i;
                int const vj = j + e1 * offset_j;
                
                std::complex<TT> const w(ptwid1[u].real(), msign * ptwid1[u].imag());
                vdst(vj) = w * vsrc(vi);
            }
            i++;
            j++;
        }
        j += jinc;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform_power_of_2(VecView<std::complex<TT>> vv)
{
    // Special case for power of 2 data
    // Process in place
    
    int logn = core_integer_log2(m_n);
    int const sign = (int)m_direction;
    
    // First use bit reversal to swap data order
    int j = 0;
    for(int i = 0; i < m_n - 1; i++)
    {
        int k = m_n / 2;
        
        if(i < j)
        {
            std::swap(vv(i), vv(j));
        }
        
        while(k <= j)
        {
            j = j - k;
            k = k / 2;
        }
        
        j += k;
    }
    
    int dual = 1;
    
    for(int bit = 0; bit < logn; bit++)
    {
        std::complex<TT> w(1,0);
        
        const double theta = 2.0 * sign * CONST_PI / (2.0 * dual);
        
        TT const s = (TT)std::sin(theta);
        TT const t = (TT)std::sin(theta / 2.0);
        TT const s2 = 2.0 * t * t;
        
        for(int b = 0; b < m_n; b += 2 * dual)
        {
            int const i = b;
            int const j = b + dual;
            
            std::complex<TT> const wd = vv(j);
            vv(j) = vv(i) - wd;
            vv(i) = vv(i) + wd;
        }
        
        for(int a = 1; a < dual; a++)
        {
            // Trig recurrence (rotation by const factor)
            TT const tmp_re = w.real() - s * w.imag() - s2 * w.real();
            TT const tmp_im = w.imag() + s * w.real() - s2 * w.imag();
            w.real(tmp_re);
            w.imag(tmp_im);
            
            for(int b = 0; b < m_n; b += 2 * dual)
            {
                int const i = b + a;
                int const j = b + a + dual;

                std::complex<TT> const wd = w * vv(j);
                vv(j) = vv(i) - wd;
                vv(i) = vv(i) + wd;
            }
        }
        
        dual *= 2;
    }
}

template <typename TT> void im::FFTComplex<TT>::transform(VecView<std::complex<TT>> vv, FFTDirection direction)
{
    IM_CHECK_VALID(vv);
    
    if(m_n<2)
        return;
    
    m_direction = direction;
    
#if SPECIALCASE_POWEROF2
    if(core_is_pow2(m_n))
    {
        transform_power_of_2(vv);
        return;
    }
#endif
    
    int product = 1;
    bool data_in_user_buffer = true;
    
    for(int i=0; i<m_factors.size(); i++)
    {
        int factor = m_factors[i];
        product *= factor;
        
        VecView<std::complex<TT>> vsrc, vdst;
        
        if(data_in_user_buffer)
        {
            data_in_user_buffer = false;
            vsrc = vv;
            vdst.wrap(m_buffer);
        }
        else
        {
            data_in_user_buffer = true;
            vsrc.wrap(m_buffer);
            vdst = vv;
        }
        
        if(factor==2)
            transform_radix_2(vdst, vsrc, product, m_offsets[i]);
        else if(factor==3)
            transform_radix_3(vdst, vsrc, product, m_offsets[i]);
        else if(factor==4)
            transform_radix_4(vdst, vsrc, product, m_offsets[i]);
        else if(factor==5)
            transform_radix_5(vdst, vsrc, product, m_offsets[i]);
        else
            transform_radix_N(vdst, vsrc, product, m_offsets[i], factor);
    }
    
    if(!data_in_user_buffer)
    {
        // Copy back
        VecView<std::complex<TT>> vbuf(m_buffer);
        vv.copy_from(vbuf);
    }
}

template class im::FFTComplex<float>;
template class im::FFTComplex<double>;

template <typename TT>
void im::core_fft_forward(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> const &vsrc)
{
    vdst.copy_from(vsrc);
    core_fft_in_place_forward(vdst);
}

#define INST(TT) template void im::core_fft_forward(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> const &vsrc)
INST(float); INST(double);
#undef INST


template <typename TT>
void im::core_fft_inverse(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> const &vsrc)
{
    vdst.copy_from(vsrc);
    core_fft_in_place_inverse(vdst);
}

#define INST(TT) template void im::core_fft_inverse(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> const &vsrc)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_in_place_forward(VecView<std::complex<TT>> v)
{
    IM_CHECK_VALID(v);
    
    FFTComplex<TT> fft;
    fft.init(v.rows());
    fft.forward(v);
}

#define INST(TT) template void im::core_fft_in_place_forward(VecView<std::complex<TT>> v)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_in_place_inverse(VecView<std::complex<TT>> v)
{
    IM_CHECK_VALID(v);
    
    FFTComplex<TT> fft;
    fft.init(v.rows());
    fft.inverse(v);
    core_block_scale(v, v, std::complex<TT>(1.0/v.rows()));
}

#define INST(TT) template void im::core_fft_in_place_inverse(VecView<std::complex<TT>> v)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_forward_real(VecView<std::complex<TT>> vdst, VecView<TT> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(vsrc.rows()==2*vdst.rows());
    
    FFTComplex<TT> fft;
    fft.init(vdst.rows());
    
    int size = vdst.rows();
    
    for(int i=0; i<size; i++)
        vdst(i) = std::complex<TT>(vsrc(i*2), vsrc(i*2+1));
    
    fft.forward(vdst);
    
    double const theta = -CONST_PI / size;
    TT const sth = std::sin(0.5*theta);
    TT const dw_re = -2*sth*sth;
    TT const dw_im = std::sin(theta);
    TT w_re = 1 + dw_re;
    TT w_im = dw_im;
    
    // note that size may be odd
    int hs = size/2;
    
    for(int i=1; i<=hs; i++)
    {
        const int smi = size - i;
        
        const TT h1_re = vdst(i).real() + vdst(smi).real();
        const TT h1_im = vdst(i).imag() - vdst(smi).imag();
        const TT h2_re = vdst(i).imag() + vdst(smi).imag();
        const TT h2_im = -vdst(i).real() + vdst(smi).real();
        
        const TT u1 = w_re * h2_re;
        const TT u2 = w_im * h2_im;
        const TT u3 = w_re * h2_im + w_im * h2_re;
        
        vdst(i) = std::complex<TT>((TT)0.5 * (h1_re + u1 - u2), (TT)0.5 * (h1_im + u3));
        vdst(smi) = std::complex<TT>((TT)0.5 * (h1_re - u1 + u2), (TT)0.5 * (-h1_im + u3));
        
        // Trig recurrence
        const TT tmp = w_re;
        w_re = w_re * dw_re - w_im * dw_im + w_re;
        w_im = w_im * dw_re + tmp * dw_im + w_im;
    }
    
    vdst(0) = std::complex<TT>(vdst(0).real() + vdst(0).imag(), vdst(0).real() - vdst(0).imag());
}

#define INST(TT) template void im::core_fft_forward_real(VecView<std::complex<TT>> vdst, VecView<TT> const &vsrc)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_inverse_real(VecView<TT> vdst, VecView<std::complex<TT>> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    IM_CHECK_VALID(vdst);
    IM_CHECK_ARGS(2*vsrc.rows()==vdst.rows());

    FFTComplex<TT> fft;
    fft.init(vsrc.rows());

    int size = vsrc.rows();
    
    im::Vec<std::complex<TT>> vbuf;
    
    if(vdst.row_stride()==1)
    {
        // compatible with packed complex, so can work in place
        vbuf.wrap(size, 1, (std::complex<TT> *)vdst.ptr());
    }
    else
    {
        // need to allocate temporary
        vbuf.resize(size);
    }
    
    vbuf.copy_from(vsrc);
    
    double const theta = CONST_PI / size;
    TT const sth = std::sin(0.5*theta);
    TT const dw_re = -2*sth*sth;
    TT const dw_im = std::sin(theta);
    TT w_re = 1 + dw_re;
    TT w_im = dw_im;
    
    // note that size may be odd
    int hs = size/2;
    
    TT const scale = (TT)(0.5)/size;
    
    for(int i=1; i<=hs; i++)
    {
        int const smi = size - i;
        
        TT const h1_re = vbuf(i).real() + vbuf(smi).real();
        TT const h1_im = vbuf(i).imag() - vbuf(smi).imag();
        TT const h2_re = -vbuf(i).imag() - vbuf(smi).imag();
        TT const h2_im = vbuf(i).real() - vbuf(smi).real();
        
        TT const u1 = w_re * h2_re;
        TT const u2 = w_im * h2_im;
        TT const u3 = w_re * h2_im + w_im * h2_re;
        
        vbuf(i).real(scale * (h1_re + u1 - u2));
        vbuf(i).imag(scale * (h1_im + u3));
        
        vbuf(smi).real(scale * (h1_re - u1 + u2));
        vbuf(smi).imag(scale * (-h1_im + u3));
        
        // Trig recurrence
        TT const tmp = w_re;
        w_re = w_re * dw_re - w_im * dw_im + w_re;
        w_im = w_im * dw_re + tmp * dw_im + w_im;
    }
    
    TT const tmp = vbuf(0).real();
    vbuf(0).real(scale * (tmp + vbuf(0).imag()));
    vbuf(0).imag(scale * (tmp - vbuf(0).imag()));
    
    fft.inverse(vbuf.view());
    
    if(vdst.row_stride()!=1)
    {
        // copy out of buf
        for(int i=0; i<size; i++)
        {
            vdst(i*2) = vbuf(i).real();
            vdst(i*2+1) = vbuf(i).imag();
        }
    }
}

#define INST(TT) template void im::core_fft_inverse_real(VecView<TT> vdst, VecView<std::complex<TT>> const &vsrc)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_forward(MtxView<std::complex<TT>> mdst, MtxView<std::complex<TT>> const &msrc)
{
    mdst.copy_from(msrc);
    core_fft_in_place_forward(mdst);
}

#define INST(TT) template void im::core_fft_forward(MtxView<std::complex<TT>> mdst, MtxView<std::complex<TT>> const &msrc)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_inverse(MtxView<std::complex<TT>> mdst, MtxView<std::complex<TT>> const &msrc)
{
    mdst.copy_from(msrc);
    core_fft_in_place_inverse(mdst);
}

#define INST(TT) template void im::core_fft_inverse(MtxView<std::complex<TT>> mdst, MtxView<std::complex<TT>> const &msrc)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_in_place_forward(MtxView<std::complex<TT>> m)
{
    IM_CHECK_VALID(m);
    
    FFTComplex<TT> fft;
    
    fft.init(m.cols());
    
    for(int i=0; i<m.rows(); i++)
        fft.forward(m.row(i));
    
    fft.init(m.rows());
    
    for(int i=0; i<m.cols(); i++)
        fft.forward(m.col(i));
}

#define INST(TT) template void im::core_fft_in_place_forward(MtxView<std::complex<TT>> m)
INST(float); INST(double);
#undef INST

template <typename TT>
void im::core_fft_in_place_inverse(MtxView<std::complex<TT>> m)
{
    IM_CHECK_VALID(m);
    
    FFTComplex<TT> fft;
    
    fft.init(m.rows());
    
    for(int i=0; i<m.cols(); i++)
        fft.inverse(m.col(i));
    
    fft.init(m.cols());
    
    for(int i=0; i<m.rows(); i++)
        fft.inverse(m.row(i));
    
    core_block_scale(m, m, std::complex<TT>(1.0/(m.rows() * m.cols())));
}

#define INST(TT) template void im::core_fft_in_place_inverse(MtxView<std::complex<TT>> m)
INST(float); INST(double);
#undef INST


