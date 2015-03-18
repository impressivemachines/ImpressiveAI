//
//  solve.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_solve_h
#define Metaphor_solve_h

namespace im
{
    enum PolySolveRoots
    {
        PolySolveRootsNoRoots,
        PolySolveRootsOneRoot,
        PolySolveRootsTwoRealRoots,
        PolySolveRootsTwoComplexRoots,
        PolySolveRootsThreeRealRoots,
        PolySolveRootsOneRealRootTwoComplexRoots
    };
    
    // Compute the roots of the equation a x^2 + b x + c = 0.
    // Return vale is the number of roots which are valid.
    // mavCoefs is a length 3 vector and mavRoots is length 2
    template <typename T> PolySolveRoots core_solve_quadratic(VecView<std::complex<T> > mavRoots, const VecView<T> &mavCoefs);
    
    // Compute the roots of the equation a x^3 + b x^2 + c x + d = 0.
    // Return vale is the number of roots which are valid.
    // mavCoefs is a length 4 vector and mavRoots is length 3
    template <typename T> PolySolveRoots core_solve_cubic(VecView<std::complex<T> > mavRoots, const VecView<T> &mavCoefs);
    
    // Compute the roots of the equation a x^2 + b x + c = 0.
    // This faster routine only works for equations with two positive real roots.
    // Roots are sorted largest to smallest.
    // mavCoefs is a length 3 vector and mavRoots is length 2
    template <typename T> void core_solve_quadratic_pos_real(VecView<T> mavRoots, const VecView<T> &mavCoefs);
    
    // Compute the roots of the equation a x^3 + b x^2 + c x + d = 0.
    // This faster routine only works for equations with three positive real roots.
    // Roots are sorted largest to smallest.
    // mavCoefs is a length 4 vector and mavRoots is length 3
    template <typename T> void core_solve_cubic_pos_real(VecView<T> mavRoots, const VecView<T> &mavCoefs);

    // Find root numerically without using derivatives using Brent method
    template <typename TT> TT core_root_search(TT (*pfunc)(TT x, void *puser), void *puser, TT bracket_min, TT bracket_max);
    
    // Fits a 1D quadratic through (-1,vf(0)) (0,vf(1)), (1,vf(2))
    // The quadratic is f(x) = coef(0) * x*x + coef(1) * x + coef(2)
    template <typename TT> void core_quadratic_fit_1d(VecView<TT> vCoefs, VecView<TT> const &vf);
    
    // Fits a 2D quadratic surface through the 3x3 set of values centered on (0,0) corresponding to mf(1,1). The y axis points down.
    // The quadratic equation is f(x,y) = 0.5 * (x,y)' (a b; b c) (x,y) + (d e)' (x,y) + f
    // Here coefs(0) through coefs(5) correspond to a through f.
    template <typename TT> void core_quadratic_fit_2d(VecView<TT> vCoefs, MtxView<TT> const &mf);
    
    // Fits a 1D quadratic and returns the position of the minimum or maximum and the function value at that point
    // The zero position corresponds to the middle sample of the length 3 vector vf
    // The function returns 1 for a maximum, -1 for a minimum, or 0 for any other condition.
    template <typename TT> int core_quadratic_interp_1d(TT &xminmax, TT &fxminmax, VecView<TT> const &vf);
    
    // Fits a 2D quadratic and returns the position of the minimum or maximum and the function value at that point
    // The zero position corresponds to the middle sample of the 3x3 matrix mf
    // The function returns 1 for a maximum, -1 for a minimum, or 0 for any other condition, e.g. flat, negligable curvature, or saddle point.
    template <typename TT> int core_quadratic_interp_2d(TT &xminmax, TT &yminmax, TT &fxminmax, MtxView<TT> const &mf);
    
    // Evaluates a 1D quadratic at the point x.
    // The quadratic is f(x) = coef(0) * x*x + coef(1) * x + coef(2)
    template <typename TT> TT core_quadratic_sample_1d(VecView<TT> const &vCoefs, TT x)
    {
        IM_CHECK_VALID(vCoefs);
        IM_CHECK_VECTOR_SIZE(vCoefs, 3);
        return x * (vCoefs(0) * x + vCoefs(1)) + vCoefs(2);
    }
    
    // Evaluates a 2D quadratic at the point x,y.
    // The quadratic equation is f(x,y) = 0.5 * (x,y)' (a b; b c) (x,y) + (d e)' (x,y) + f
    // Here coefs(0) through coefs(5) correspond to a through f.
    template <typename TT> TT core_quadratic_sample_2d(VecView<TT> const &vCoefs, TT x, TT y)
    {
        IM_CHECK_VALID(vCoefs);
        IM_CHECK_VECTOR_SIZE(vCoefs, 6);
        return (TT)0.5 * (vCoefs(0) * x*x + vCoefs(1) * x*y + vCoefs(2) * y*y) + vCoefs(3) * x + vCoefs(4) * y + vCoefs(5);
    }
}

#endif
