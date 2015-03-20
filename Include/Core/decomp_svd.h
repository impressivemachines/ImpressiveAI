//
//  decomp_svd.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_decomp_svd_h
#define Metaphor_decomp_svd_h

namespace im
{
    // Decompose an arbitrary m x n matrix using SVD:
    // A = USV'
    // Where U is an m x n matrix, V is n x n unitary matrix, and S is n x n real positive matrix with zero outside the diagonal
    // The singular values are the diagonal entries of S and are sorted into decreasing order.
    
    // This class is only implemented for m >= n.
    
    template <typename TT>
    class MatrixDecompSVD
    {
    public:
        MatrixDecompSVD() {}
        MatrixDecompSVD(MtxView<TT> const &mA) { compute(mA); }
        
        void compute(MtxView<TT> const &mA, bool qr_precondition = true);
        
        // When a singular value sv is <= sv_thresh then 1/sv is set to zero when forming the solution
        // This deals with near singularity of A
        Mtx<TT> solve(MtxView<TT> const &mavy, TT sv_thresh = (TT)0) const;
        Vec<TT> solve(VecView<TT> const &vvy, TT sv_thresh = (TT)0) const;
        
        Mtx<TT> pseudo_inverse(TT sv_thresh = (TT)0) const;
        
        Vec<TT> const vectorS() const { return m_vS; }
        Mtx<TT> const matrixU() const { return m_mU; }
        Mtx<TT> const matrixV() const { return m_mV; }
        
        void deallocate() { m_vS.deallocate(); m_mU.deallocate(); m_mV.deallocate(); }
        
        // Compute SVD in place without preconditioning
        void in_place_jacobi(MtxView<TT> mvA, MtxView<TT> mvV, VecView<TT> vvS) const;
   
    private:
        Vec<TT> m_vS;
        Mtx<TT> m_mU;
        Mtx<TT> m_mV;
    };
    
    // Special case of 2x2 SVD
    template <typename TT> void core_decomp_svd_2x2(MtxView<TT> mvU, VecView<TT> vvS, MtxView<TT> mvV, MtxView<TT> const &mvA);
    
    
}

#endif
