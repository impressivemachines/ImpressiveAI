//
//  solve.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

// Compute the roots of the equation a x^2 + b x + c = 0.
// Return vale is the number of roots which are valid.
template <typename T> im::PolySolveRoots im::core_solve_quadratic(VecView<std::complex<T> > vRoots, const VecView<T> &vCoefs)
{
    IM_CHECK_VALID(vRoots);
    IM_CHECK_VALID(vCoefs);
    IM_CHECK_VECTOR_SIZE(vRoots, 2);
    IM_CHECK_VECTOR_SIZE(vCoefs, 3);
    
    T a = vCoefs(0);
    T b = vCoefs(1);
    T c = vCoefs(2);
    
    if(a==0 && b==0)
    {
        vRoots(0) = 0;
        vRoots(1) = 0;
        
        return im::PolySolveRootsNoRoots;
    }
    
    if(a==0)
    {
        vRoots(0) = std::complex<T>(-c/b);
        vRoots(1) = 0;
        
        return im::PolySolveRootsOneRoot;
    }
    
    T s = b*b - 4*a*c;
    
    if(s<0)
    {
        // complex vRoots
        T re = -b/(2*a);
        T im = std::sqrt(-s)/(2*a);
        vRoots(0) = std::complex<T>(re, im);
        vRoots(1) = std::complex<T>(re, -im);
        
        return im::PolySolveRootsTwoComplexRoots;
    }
    
    // real roots
    T q;
    if(b>=0)
        q = -(T)0.5 * (b + std::sqrt(s));
    else
        q = -(T)0.5 * (b - std::sqrt(s));
    
    vRoots(0) = std::complex<T>(q/a);
    vRoots(1) = std::complex<T>(c/q);
    
    return im::PolySolveRootsTwoRealRoots;
}

template im::PolySolveRoots im::core_solve_quadratic(VecView<im::Cf > vRoots, const VecView<float> &vCoefs);
template im::PolySolveRoots im::core_solve_quadratic(VecView<im::Cd > vRoots, const VecView<double> &vCoefs);

// Compute the roots of the equation a x^3 + b x^2 + c x + d = 0.
// Return vale is the number of roots which are valid.
template <typename T> im::PolySolveRoots im::core_solve_cubic(VecView<std::complex<T> > vRoots, const VecView<T> &vCoefs)
{
    IM_CHECK_VALID(vRoots);
    IM_CHECK_VALID(vCoefs);
    IM_CHECK_VECTOR_SIZE(vRoots, 3);
    IM_CHECK_VECTOR_SIZE(vCoefs, 4);
    
    T s = vCoefs(0);
    if(s==0)
    {
        im::PolySolveRoots result = core_solve_quadratic<T>(vRoots.head(2), vCoefs.tail(3));
        vRoots(2) = 0;
        
        return result;
    }
    
    T a = vCoefs(1);
    T b = vCoefs(2);
    T c = vCoefs(3);
    
    if(s!=1)
    {
        a/=s;
        b/=s;
        c/=s;
    }
    
    T aa = a*a;
    T q = (aa - 3*b)/9;
    T r = (2*a*aa - 9*a*b + 27*c)/54;
    T rr = r*r;
    T qqq = q*q*q;
    
    if(rr < qqq)
    {
        T th = acos(r/std::sqrt(qqq));
        T mag = -2*std::sqrt(q);
        T a3 = a/3;
        
        vRoots(0) = std::complex<T>(mag*std::cos(th/3) - a3);
        vRoots(1) = std::complex<T>(mag*std::cos((th + 2*CONST_PI)/3) - a3);
        vRoots(2) = std::complex<T>(mag*std::cos((th - 2*CONST_PI)/3) - a3);
        
        return im::PolySolveRootsThreeRealRoots;
    }
    
    T w = std::pow(fabs(r) + std::sqrt(rr - qqq), 1/(T)3.0);
    if(w * r > 0)
        w *= -1;
    T f = 0;
    if(w!=0)
        f = q/w;
    T u = -(T)0.5*(w + f) - a/3;
    T v = (T)(CONST_SQRT3_2) * (w - f);
    
    vRoots(0) = std::complex<T>(w + f - a/(T)3.0);
    vRoots(1) = std::complex<T>(u, v);
    vRoots(2) = std::complex<T>(u, -v);
    
    return im::PolySolveRootsOneRealRootTwoComplexRoots;
}

template im::PolySolveRoots im::core_solve_cubic(VecView<im::Cf> vRoots, const VecView<float> &vCoefs);
template im::PolySolveRoots im::core_solve_cubic(VecView<im::Cd> vRoots, const VecView<double> &vCoefs);

// Compute the roots of the equation a x^2 + b x + c = 0.
// This routine only works for equations with two positive real roots.
template <typename T> void im::core_solve_quadratic_pos_real(VecView<T> vRoots, const VecView<T> &vCoefs)
{
    IM_CHECK_VALID(vRoots);
    IM_CHECK_VALID(vCoefs);
    IM_CHECK_VECTOR_SIZE(vRoots, 2);
    IM_CHECK_VECTOR_SIZE(vCoefs, 3);
    
    T b = -vCoefs(1) / vCoefs(0);
    T c = vCoefs(2) / vCoefs(0);
    
    T b2 = (T)0.5 * b;
    T d = b2 * b2 - c;
    d = (d > 0) ? (T)std::sqrt(d) : (T)0.0;
    
    vRoots(0) = std::min(b2 + d, b);
    vRoots(1) = std::max(b2 - d, (T)0);
}

template void im::core_solve_quadratic_pos_real(VecView<float> vRoots, const VecView<float> &vCoefs);
template void im::core_solve_quadratic_pos_real(VecView<double> vRoots, const VecView<double> &vCoefs);

// Compute the roots of the equation a x^3 + b x^2 + c x + d = 0.
// This routine only works for equations with three positive real roots.
template <typename T> void im::core_solve_cubic_pos_real(VecView<T> vRoots, const VecView<T> &vCoefs)
{
    IM_CHECK_VALID(vRoots);
    IM_CHECK_VALID(vCoefs);
    IM_CHECK_VECTOR_SIZE(vRoots, 3);
    IM_CHECK_VECTOR_SIZE(vCoefs, 4);
    
    T s = vCoefs(0);
    T a = vCoefs(1);
    T b = vCoefs(2);
    T c = vCoefs(3);
    
    if(s!=1)
    {
        a/=s;
        b/=s;
        c/=s;
    }
    
    T aa = a*a;
    T q = (aa - 3*b)/9;
    T r = (2*a*aa - 9*a*b + 27*c)/54;
    T rr = r*r;
    T qqq = q*q*q;
    
    if(rr < qqq)
    {
        T a3 = a/3;
        vRoots(0) = -a3;
        vRoots(1) = -a3;
        vRoots(2) = -a3;
        
        if(q<=0)
            return;
        
        T mag = 2*std::sqrt(q);
        T th = acos(r/std::sqrt(qqq));
        
        vRoots(0) -= mag * std::cos(th/3);
        vRoots(1) -= mag * std::cos((th + 2*CONST_PI)/3);
        vRoots(2) -= mag * std::cos((th - 2*CONST_PI)/3);
        
        // sort largest to smallest
        core_sort_3(vRoots(0), vRoots(1), vRoots(2), SortDirectionDescending);
    }
    else
    {
        vRoots(0) = 0;
        vRoots(1) = 0;
        vRoots(2) = 0;
    }
}

template void im::core_solve_cubic_pos_real(VecView<float> vRoots, const VecView<float> &vCoefs);
template void im::core_solve_cubic_pos_real(VecView<double> vRoots, const VecView<double> &vCoefs);

template <typename TT>
TT im::core_root_search(TT (*pfunc)(TT x, void *puser), void *puser, TT bracket_min, TT bracket_max)
{
    IM_CHECK_ARGS(bracket_max >= bracket_min);
    
    if(bracket_max==bracket_min)
        return bracket_min;
    
    TT a = bracket_min;
    TT fa = (*pfunc)(a, puser);
    TT b = bracket_max;
    TT fb = (*pfunc)(b, puser);
    
    if((fa<(TT)0 && fb<(TT)0) || (fa>(TT)0 && fb>(TT)0))
        IM_THROW_NO_SOLUTION; // root is not in this interval
    
    TT c = a;
    TT fc = fa;
    TT d = b-a;
    TT e = b-a;
    
    for(int i=0; i<1000; i++)
    {
        if(std::abs(fc) < std::abs(fb))
        {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        
        TT tol = (TT)0.5 * TypeProperties<TT>::epsilon() * std::abs(b);
        TT m = (TT)0.5 * (c-b);
        
        if(fb==(TT)0 || std::abs(m) <= tol)
            return b;
        
        if(std::abs(e) < tol || std::abs(fa) <= std::abs(fb))
        {
            // use bisection
            d = m;
            e = m;
        }
        else
        {
            // try interpolation
            TT s = fb / fa;
            TT p,q,r;
            
            if(a==c)
            {
                p = (TT)2 * m * s;
                q = (TT)1 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * ((TT)2 * m * q * (q-r) - (b-a) * (r-(TT)1));
                q = (q-(TT)1) * (r-(TT)1) * (s-(TT)1);
            }
            
            if(p>(TT)0)
                q = -q;
            else
                p = -p;
            
            TT p2 = (TT)2 * p;
            if(p2 < std::abs(e * q) && p2 < (TT)3 * m * q - std::abs(tol * q))
            {
                // interpolation ok
                e = d;
                d = p/q;
            }
            else
            {
                // fall back on bisection
                d = m;
                e = m;
            }
        }
        
        a = b;
        fa = fb;
        
        if(tol < std::abs(d))
            b += d;
        else
            b += (m>(TT)0 ? tol : -tol);
        
        fb = (*pfunc)(b, puser);
        
        if((fb < (TT)0 && fc < (TT)0) || (fb > (TT)0 && fc > (TT)0))
        {
            c = a;
            fc = fa;
            d = e = b-a;
        }
    }
    
    IM_THROW_NO_SOLUTION;
    
    return (TT)0;
}

template float im::core_root_search(float (*pfunc)(float x, void *puser), void *puser, float bracket_min, float bracket_max);
template double im::core_root_search(double (*pfunc)(double x, void *puser), void *puser, double bracket_min, double bracket_max);

template <typename TT> void im::core_quadratic_fit_1d(VecView<TT> vCoefs, VecView<TT> const &vf)
{
    IM_CHECK_VALID(vf);
    IM_CHECK_VALID(vCoefs);
    IM_CHECK_VECTOR_SIZE(vf, 3);
    IM_CHECK_VECTOR_SIZE(vCoefs, 3);
    
    vCoefs(0) = (TT)0.5 * (vf(2) + vf(0)) - vf(1);
    vCoefs(1) = (TT)0.5 * (vf(2) - vf(0));
    vCoefs(2) = vf(0);
}

template void im::core_quadratic_fit_1d(VecView<float> vCoefs, VecView<float> const &vf);
template void im::core_quadratic_fit_1d(VecView<double> vCoefs, VecView<double> const &vf);

template <typename TT> void im::core_quadratic_fit_2d(VecView<TT> vCoefs, MtxView<TT> const &mf)
{
    IM_CHECK_VALID(mf);
    IM_CHECK_VALID(vCoefs);
    IM_CHECK_MATRIX_SIZE(mf, 3, 3);
    IM_CHECK_VECTOR_SIZE(vCoefs, 6);
    
    vCoefs(0) = (mf(0,0) + mf(1,0) + mf(2,0) - (TT)2 * (mf(0,1) + mf(1,1) + mf(2,1)) + mf(0,2) + mf(1,2) + mf(2,2)) / (TT)3;
    vCoefs(1) = (mf(0,0) - mf(0,2) - mf(2,0) + mf(2,2))/(TT)4;
    vCoefs(2) = (mf(0,0) + mf(0,1) + mf(0,2) - (TT)2 * (mf(1,0) + mf(1,1) + mf(1,2)) + mf(2,0) + mf(2,1) + mf(2,2)) / (TT)3;
    vCoefs(3) = (mf(0,2) + mf(1,2) + mf(2,2) - mf(0,0) - mf(1,0) - mf(2,0))/(TT)6;
    vCoefs(4) = (mf(2,0) + mf(2,1) + mf(2,2) - mf(0,0) - mf(0,1) - mf(0,2))/(TT)6;
    vCoefs(5) = ((TT)2 * (mf(0,1) + mf(1,0) + mf(1,2) + mf(2,1)) + (TT)5 * mf(1,1) - mf(0,0) - mf(0,2) - mf(2,0) - mf(2,2))/(TT)9;
}

template void im::core_quadratic_fit_2d(VecView<float> vCoefs, MtxView<float> const &mf);
template void im::core_quadratic_fit_2d(VecView<double> vCoefs, MtxView<double> const &mf);

template <typename TT> int im::core_quadratic_minmax_1d(TT &xminmax, TT &fxminmax, VecView<TT> const &vf)
{
    TT coefs[3];
    core_quadratic_fit_1d(VecView<TT>(3, 1, coefs), vf);
 
    if(std::abs(coefs[0]) < TypeProperties<TT>::epsilon())
    {
        xminmax = (TT)0;
        fxminmax = coefs[2];
        return 0;
    }
    
    xminmax = -coefs[1]/((TT)2*coefs[0]);
    fxminmax = core_quadratic_sample_1d(VecView<TT>(3, 1, coefs), xminmax);
    return coefs[0] > (TT)0 ? -1 : 1;
}

template int im::core_quadratic_minmax_1d(float &xminmax, float &fxminmax, VecView<float> const &vf);
template int im::core_quadratic_minmax_1d(double &xminmax, double &fxminmax, VecView<double> const &vf);

template <typename TT> int im::core_quadratic_minmax_2d(TT &xminmax, TT &yminmax, TT &fxminmax, MtxView<TT> const &mf)
{
    TT coefs[6];
    core_quadratic_fit_2d(VecView<TT>(6, 1, coefs), mf);
    
    // determinant is product of eigenvalues, so if negative they have opposite signs
    TT det = coefs[0] * coefs[2] - coefs[1] * coefs[1];
    if(det <= TypeProperties<TT>::epsilon())
    {
        // not sufficiently curved or else a saddle point
        xminmax = (TT)0;
        yminmax = (TT)0;
        fxminmax = coefs[5];
        return 0;
    }
    
    xminmax = (-coefs[2] * coefs[3] + coefs[1] * coefs[4]) / det;
    yminmax = (coefs[1] * coefs[3] - coefs[0] * coefs[4]) / det;
    fxminmax = core_quadratic_sample_2d(VecView<TT>(6, 1, coefs), xminmax, yminmax);

    return coefs[0] > (TT)0 ? -1 : 1;
}

template int im::core_quadratic_minmax_2d(float &xminmax, float &yminmax, float &fxminmax, MtxView<float> const &mf);
template int im::core_quadratic_minmax_2d(double &xminmax, double &yminmax, double &fxminmax, MtxView<double> const &mf);

