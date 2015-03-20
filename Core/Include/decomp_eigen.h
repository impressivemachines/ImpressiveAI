//
//  decomp_eigen.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_decomp_eigen_h
#define Metaphor_decomp_eigen_h

namespace im
{
    // Computes the eigen decomposition of a symmetric matrix A according to
    // A = V D V' where D is diagonal and V is orthogonal with the eigenvectors as its columns.
    // The diagonal elements of D are the eigenvalues.
    // Since A is symmetric, the eigenvalues are all positive and real.
    
    template <typename TT>
    class MatrixDecompEigenSymmetric
    {
    public:
        MatrixDecompEigenSymmetric() : m_vectors(false) {}
        MatrixDecompEigenSymmetric(MtxView<TT> const &mavA) : m_vectors(false) { compute(mavA); }
        
        void compute(MtxView<TT> const &mavA, bool compute_vectors = true);
        
        Mtx<TT> solve(MtxView<TT> const &mavy) const;
        Vec<TT> solve(VecView<TT> const &vvy) const;
        
        // Get the 2D matrix V with eigenvectors as the columns
        Mtx<TT> const eigenvectors() const { return m_mvectors; }
        
        // Get the 1D vector of eigenvalues corresponding to the columns of V
        Vec<TT> const eigenvalues() const { return m_vdiag; }
        
        // Computes V D^(1/2) V'
        Mtx<TT> square_root() const;
        
        // Computes V D^(-1/2) V'
        Mtx<TT> inverse_square_root() const;
        
        void deallocate() { m_mvectors.deallocate(); m_vdiag.deallocate(); m_vsubdiag.deallocate(); }
        
    private:
        Mtx<TT> m_mvectors;
        Vec<TT> m_vdiag;
        Vec<TT> m_vsubdiag;
        bool m_vectors;
    };
    
    // Special case 2x2 symmetric eigen decomposition
    template <typename TT> void core_decomp_eigen_2x2(VecView<TT> vevals, MtxView<TT> mevecs, MtxView<TT> const &mA);
    
    
    
    // TODO: general eigen decomp
    
}

#endif
