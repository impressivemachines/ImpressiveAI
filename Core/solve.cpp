//
//  solve.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

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

