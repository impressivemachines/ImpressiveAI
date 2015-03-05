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

    
}

#endif
