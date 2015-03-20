//
//  decomp_tools.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_decomp_tools_h
#define Metaphor_decomp_tools_h

namespace im
{
    // Transforms the given column vector in place into a householder vector and returns the coefficient that
    // would annihilate elements 1 through n-1
    template <typename TT> TT core_compute_householder_vector(VecView<TT> vhh);
    
    // Applies a householder transform vector and tau to the matrix
    // A = A - tau w v' to annihilate column coefficients below the diagonal
    template <typename TT> void core_apply_householder_vector(MtxView<TT> mavA, TT const &tau, VecView<TT> const &vhh);
    
    // Applies a householder transformation to a matrix being build up from the identity matrix, using the first column of A as
    // a householder vector
    template <typename TT> void core_apply_householder_identity(MtxView<TT> mavA, TT const &tau);
}

#endif
