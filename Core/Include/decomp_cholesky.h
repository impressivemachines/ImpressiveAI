//
//  decomp_cholesky.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_decomp_cholesky_h
#define Metaphor_decomp_cholesky_h

namespace im
{
    // Decomposition of **positive definite symmetric** matrices
    
    // LDL decomposition where A = L D L^T, L is lower triangular with 1s on the diagonal, and D is diagonal.
    
    // TODO: support complex matrices
    
    template <typename TT>
    class MatrixDecompLDLT
    {
    public:
        MatrixDecompLDLT() {}
        MatrixDecompLDLT(MtxView<TT> const &mavA) { compute(mavA); } // A must be positive definite symmtric
        
        // only the lower triangular section is referenced
        void compute(MtxView<TT> const &mavA);
        
        // Solve A X = Y
        Mtx<TT> solve(MtxView<TT> const &mavy) const;
        Vec<TT> solve(VecView<TT> const &mavy) const;
        
        // Get the lower triangular matrix
        Mtx<TT> matrix_L() const
        {
            Mtx<TT> matrtn(m_matLDLT.rows(), m_matLDLT.cols());
            core_block_clear_upper_tri(matrtn.view(), false);
            matrtn.diag() = (TT)1;
            core_block_copy_lower_tri(matrtn.view(), m_matLDLT.view(), false);
            return matrtn;
        }
        
        // Get a 1D vector for the diagonal elements
        Vec<TT> const vector_D() const { return m_matLDLT.diag(); }
        
        // Compute the inverse
        Mtx<TT> inverse() const;
        
        // Get the symmetrical packed LDLT matrix
        Mtx<TT> const matrix_LDLT() const { return m_matLDLT; }
        
        void deallocate() { m_matLDLT.deallocate(); }
        
    private:
        Mtx<TT> m_matLDLT;
    };
    
    // Cholesky decomposition where A = L L^T, L is lower triangular.
    // A is positive semidef symmetric
    // Matrix square root
    
    template <typename TT>
    class MatrixDecompLLT
    {
    public:
        MatrixDecompLLT() {}
        MatrixDecompLLT(MtxView<TT> const &mavA) { compute(mavA); } // A must be positive definite symmtric
        
        // only the lower triangular section is referenced
        void compute(MtxView<TT> const &mavA);
        
        // Solve A X = Y
        // X and Y can have multiple columns to solve for multiple vectors
        Mtx<TT> solve(MtxView<TT> const &mavy) const;
        Vec<TT> solve(VecView<TT> const &mavy) const;
        
        // Get the lower triangular matrix
        Mtx<TT> matrix_L() const
        {
            Mtx<TT> matrtn(m_matLLT.rows(), m_matLLT.cols());
            core_block_clear_upper_tri(matrtn.view(), false);
            core_block_copy_lower_tri(matrtn.view(), m_matLLT.view(), true);
            return matrtn;
        }
        
        // Compute the inverse
        Mtx<TT> inverse() const;
        
        // Get the symmetrical packed LLT matrix
        Mtx<TT> const matrix_LLT() const { return m_matLLT; }

    private:
        Mtx<TT> m_matLLT;
    };
}

#endif
