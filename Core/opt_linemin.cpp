//
//  opt_linemin.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/10/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"


template <typename TT>
void im::core_line_min_bracket(TT &xa, TT &xb, TT &xc, FuncEval1D<TT> *peval, TT bmin, TT bmax)
{
    TT const golden = (TT)CONST_G;
    TT const eps = TypeProperties<TT>::epsilon();
    
    xa = bmin;
    TT fa = peval->eval_fx(xa);
    xb = bmax;
    TT fb = peval->eval_fx(xb);
    
    if(fb>fa)
    {
        // maintain fa as largest
        std::swap(xa, xb);
        std::swap(fa, fb);
    }
    
    xc = xb + golden * (xb-xa); // project downhill
    TT fc = peval->eval_fx(xc);

    if(fb<fa && fb<fc)
    {
        if(xa>xc)
            std::swap(xa, xc);
        return;
    }
    
    int count;
    int const countmax = 100;
    
    for(count = 0; count<countmax && fb > fc; count++)
    {
        // try parabolic fit to extrapolate using abc
        TT p1 = (xb-xa) * (fb-fc);
        TT p2 = (xb-xc) * (fb-fa);
        TT denom = (TT)2 * std::max(std::abs(p2-p1), eps);
        if(p2<p1)
            denom = -denom;
        
        // new trial point
        TT xx = xb - ((xb-xc)*p2 - (xb-xa)*p1) / denom;
        
        // upper limit we want to risk
        TT xmax = xb + (TT)100 * (xc-xb);
        
        TT fx;
        
        if((xb-xx)*(xx-xc) > (TT)0)
        {
            // x lies between b and c
            fx = peval->eval_fx(xx);
            if(fx < fc)
            {
                // min between b and c
                xa = xb;
                xb = xx;
                if(xa>xc)
                    std::swap(xa, xc);
                return;
            }
            if(fx > fb)
            {
                // min between a and x
                xc = xx;
                if(xa>xc)
                    std::swap(xa, xc);
                return;
            }
            
            // fall back on golden section
            xx = xc + golden * (xc-xb);
            fx = peval->eval_fx(xx);
        }
        else if((xc-xx)*(xx-xmax) > (TT)0)
        {
            // x lies between c and xmax
            fx = peval->eval_fx(xx);
            if(fx < fc)
            {
                TT tmp = xx + golden * (xx-xc);
                xb = xc;
                xc = xx;
                xx = tmp;
                fb = fc;
                fc = fx;
                fx = peval->eval_fx(xx);
            }
        }
        else if((xx-xmax)*(xmax-xc) >= (TT)0)
        {
            xx = xmax;
            fx = peval->eval_fx(xx);
        }
        else
        {
            // fall back on golden section
            xx = xc + golden * (xc-xb);
            fx = peval->eval_fx(xx);
        }

        xa = xb;
        xb = xc;
        xc = xx;
        fa = fb;
        fb = fc;
        fc = fx;
    }
    
    if(xa>xc)
        std::swap(xa, xc);
    
    if(count==countmax)
        IM_THROW_NO_SOLUTION;

    return;
}

template void im::core_line_min_bracket(float &xa, float &xb, float &xc, FuncEval1D<float> *peval, float bmin, float bmax);
template void im::core_line_min_bracket(double &xa, double &xb, double &xc, FuncEval1D<double> *peval, double bmin, double bmax);

template <typename TT>
void im::core_line_min(TT &xmin, TT &fxmin, FuncEval1D<TT> *peval, TT bmin, TT bmax)
{
    TT const golden = (TT)(1.0  / (CONST_G + 1.0));

    if(bmin>bmax)
        std::swap(bmin, bmax);
    
    TT xa = bmin;
    TT xb = bmax;
    
    TT x = xa + golden * (xb - xa);
    TT w = x;
    TT v = x;
    
    TT fx = peval->eval_fx(x);
    TT fw = fx;
    TT fv = fx;
    
    TT d = (TT)0;
    TT e = (TT)0;
    
    TT const eps = TypeProperties<TT>::epsilon();
    TT const t = TypeProperties<TT>::epsilon() * (TT)1.0e-3;
    
    int count;
    for(count = 0; count<100; count++)
    {
        TT m = (TT)0.5 * (xa + xb);
        TT tol1 = eps * std::abs(x) + t;
        TT tol2 = (TT)2 * tol1;
        
        if(std::abs(x-m) <= tol2 - (TT)0.5 * (xb-xa))
        {
            xmin = x;
            fxmin = fx;
            return;
        }
        
        bool going_for_gold = true;
        
        if(std::abs(e) > tol1)
        {
            // fit parabola
            TT r = (x-w) * (fx-fv);
            TT q = (x-v) * (fx-fw);
            TT p = (x-v) * q - (x-w) * r;
            q = (TT)2 * (q-r);
            
            if(q>(TT)0)
                p = -p;
            else
                q = -q;
            
            r = e;
            e = d;
            
            if(std::abs(p) < std::abs((TT)0.5 * q * r) && p > q*(xa-x) && p < q*(xb-x))
            {
                // use the parabolic step
                going_for_gold = false;
                d = p / q;
                TT u = x + d;
                if((u-xa)<tol2 || (xb-u)<tol2)
                {
                    if(x<m)
                        d = tol1;
                    else
                        d = -tol1;
                }
            }
        }
        
        if(going_for_gold)
        {
            // default to golden section
            if(x<m)
                e = xb - x;
            else
                e = xa - x;
            d = golden * e;
        }
        
        TT u;
        if(std::abs(d) >= tol1)
            u = x + d;
        else if(d>(TT)0)
            u = x + tol1;
        else
            u = x - tol1;
        
        TT fu = peval->eval_fx(u);
        
        if(fu<=fx)
        {
            if(u<x)
                xb = x;
            else
                xa = x;
            
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if(u<x)
                xa = u;
            else
                xb = u;
            if(fu<=fw || w==x)
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if(fu<=fv || v==x || v==w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    
    xmin = x;
    fxmin = fx;

    IM_THROW_NO_SOLUTION;
}

template void im::core_line_min(float &xmin, float &fxmin, FuncEval1D<float> *peval, float bmin, float bmax);
template void im::core_line_min(double &xmin, double &fxmin, FuncEval1D<double> *peval, double bmin, double bmax);

template <typename TT>
void im::core_line_min_using_derivs(TT &xmin, TT &fxmin, FuncEval1D<TT> *peval, TT bmin, TT bmax)
{
    if(bmin>bmax)
        std::swap(bmin, bmax);
    
    TT xa = bmin;
    TT xb = bmax;
    
    TT x = (TT)0.5 * (xa + xb);
    TT w = x;
    TT v = x;
    
    TT fx = peval->eval_fx(x);
    TT fw = fx;
    TT fv = fx;
    
    TT dfx = peval->eval_dfx(x);
    TT dfw = dfx;
    TT dfv = dfx;
    
    TT d = (TT)0;
    TT e = (TT)0;
    
    TT const eps = TypeProperties<TT>::epsilon();
    TT const t = TypeProperties<TT>::epsilon() * (TT)1.0e-3;
    
    int count;
    for(count = 0; count<100; count++)
    {
        TT m = (TT)0.5 * (xa + xb);
        TT tol1 = eps * std::abs(x) + t;
        TT tol2 = (TT)2 * tol1;
        
        if(std::abs(x-m) <= tol2 - (TT)0.5 * (xb-xa))
        {
            xmin = x;
            fxmin = fx;
            return;
        }

        bool going_for_bisect = true;

        if(std::abs(e) > tol1)
        {
            TT d1 = (TT)2 * (xb-xa);
            TT d2 = d1;
            
            // secant method
            if(dfw != dfx)
                d1 = (w-x) * dfx / (dfx-dfw);
            if(dfv != dfx)
                d2 = (v-x) * dfx / (dfx-dfv);
            
            TT u1 = x + d1;
            TT u2 = x + d2;
            
            bool good1 = (xa-u1)*(u1-xb) > (TT)0 && dfx * d1 <= (TT)0;
            bool good2 = (xa-u2)*(u2-xb) > (TT)0 && dfx * d2 <= (TT)0;
            
            TT r = e;
            e = d;
            
            if(good1 || good2)
            {
                if(good1 && good2)
                {
                    if(std::abs(d1) < std::abs(d2))
                        d = d1;
                    else
                        d = d2;
                }
                else if(good1)
                    d = d1;
                else
                    d = d2;
                
                if(std::abs(d) <= std::abs((TT)0.5 * r))
                {
                    // ok to use secant
                    going_for_bisect = false;
                    
                    TT u = x + d;
                    if((u-xa)<tol2 || (xb-u)<tol2)
                    {
                        if(x<m)
                            d = tol1;
                        else
                            d = -tol1;
                    }
                }
            }
        }
        
        if(going_for_bisect)
        {
            // have to default to bisection
            if(dfx>=(TT)0)
                e = xa - x;
            else
                e = xb - x;
            
            d = (TT)0.5 * e;
        }
    
        TT u, fu;
        if(std::abs(d) >= tol1)
        {
            u = x + d;
            fu = peval->eval_fx(u);
        }
        else
        {
            if(d>(TT)0)
                u = x + tol1;
            else
                u = x - tol1;
            
            fu = peval->eval_fx(u);
            
            if(fu > fx)
            {
                xmin = x;
                fxmin = fx;
                return;
            }
        }
        
        TT dfu = peval->eval_dfx(u);
        
        if(fu<=fx)
        {
            if(u<x)
                xb = x;
            else
                xa = x;
            
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
            dfv = dfw;
            dfw = dfx;
            dfx = dfu;
        }
        else
        {
            if(u<x)
                xa = u;
            else
                xb = u;
            if(fu<=fw || w==x)
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
                dfv = dfw;
                dfw = dfu;
            }
            else if(fu<=fv || v==x || v==w)
            {
                v = u;
                fv = fu;
                dfv = dfu;
            }
        }
    }
    
    xmin = x;
    fxmin = fx;
    
    IM_THROW_NO_SOLUTION;
}

template void im::core_line_min_using_derivs(float &xmin, float &fxmin, FuncEval1D<float> *peval, float bmin, float bmax);
template void im::core_line_min_using_derivs(double &xmin, double &fxmin, FuncEval1D<double> *peval, double bmin, double bmax);
