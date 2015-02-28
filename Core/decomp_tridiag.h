//
//  decomp_tridiag.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_decomp_tridiag_h
#define ImpressiveAI_decomp_tridiag_h

namespace im
{
    // Reduction of a square self-adjoint matrix to tri-diagonal form
    // A = Q T Q^T, where Q is orthogonal and T is symmetric tri-diagonal
    // Only the diagonal and lower triangular part of A is referenced
    
    template <typename TT>
    class MatrixDecompTridiag
    {
    public:
        MatrixDecompTridiag() {}
        MatrixDecompTridiag(MtxView<TT> const &mavA) { compute(mavA); }
        
        void compute(MtxView<TT> const &mavA);
        
        void compute_in_place(MtxView<TT> mavA, VecView<TT> vTau); // for pro users:)
        
        // Gets the matrix Q
        Mtx<TT> matrix_Q() const;
        
        // Gets the matrix T
        Mtx<TT> matrix_T() const;
        
        // Gets the diagonal of the matrix T
        Vec<TT> const diagonal() const { return m_hh.diag(); }
        
        // Gets the sub-diagonal of the matrix T
        Vec<TT> const sub_diagonal() const { return m_hh.diag(-1); }
        
        // Gets the householder coefficients
        Vec<TT> const householder_coefs() const { return m_vtau; }
        
        // Gets the packed form of the householder matrix
        Mtx<TT> const householder_packed() const { return m_hh; }
        
        void deallocate() { m_vtau.deallocate(); m_hh.deallocate(); }
        
    private:
        Vec<TT> m_vtau;
        Mtx<TT> m_hh;
    };
    
}

#endif
