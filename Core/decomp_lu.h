//
//  decomp_lu.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_decomp_lu_h
#define Metaphor_decomp_lu_h

namespace im
{
    // Decomposition of general square matrix
    // P A = L U
    // where P is a permutation matrix, L is unit lower triangular with 1s on diagonal, and U is upper triangular.
    
    // Works for real and complex matrices.
    
    template <typename TT>
    class MatrixDecompLU
    {
    public:
        MatrixDecompLU() {}
        MatrixDecompLU(MtxView<TT> const &mavA) { compute(mavA); }
        
        // Compute the decomposition. A must be square
        void compute(MtxView<TT> const &mavA);
        
        // Solve AX = Y
        Mtx<TT> solve(MtxView<TT> const &mavy) const;
        Vec<TT> solve(VecView<TT> const &vvy) const;
        
        // Compute the inverse
        Mtx<TT> inverse() const;
        
        // Returns the determinant
        TT det() const
        {
            TT det = (TT)m_perm.sign();
            for(int i=0; i<m_matLU.rows(); i++)
                det *= m_matLU(i,i);
            return det;
        }
        
        // Returns the log of the abs(determinant)
        TT log_det() const
        {
            TT logdet = (TT)0;
            for(int i=0; i<m_matLU.rows(); i++)
                logdet += std::log(std::abs(m_matLU(i,i)));
            return logdet;
        }
        
        // Returns the sign of the determinant, or 0 if singular
        TT sign_det() const
        {
            TT s = (TT)m_perm.sign();
            for(int i=0; i<m_matLU.rows(); i++)
            {
                TT val = m_matLU(i,i);
                if(val==(TT)0)
                {
                    s = (TT)0;
                    break;
                }
                else
                    s *= val/std::abs(val); // because TT may be complex
            }
            return s;
        }
        
        // Test for singular
        bool is_singuar() const
        {
            for(int i=0; i<m_matLU.rows(); i++)
                if(m_matLU(i,i)==(TT)0)
                    return true;
            return false;
        }
        
        Mtx<TT> matrix_L() const
        {
            Mtx<TT> matrtn(m_matLU.rows(), m_matLU.cols());
            core_block_copy_lower_tri(matrtn.view(), m_matLU.view(), false);
            core_block_clear_upper_tri(matrtn.view(), false);
            matrtn.diag() = (TT)1;
            return matrtn;
        }
        
        Mtx<TT> matrix_U() const
        {
            Mtx<TT> matrtn(m_matLU.rows(), m_matLU.cols());
            core_block_copy_upper_tri(matrtn.view(), m_matLU.view(), true);
            core_block_clear_lower_tri(matrtn.view(), false);
            return matrtn;
        }
        
        Mtx<TT> matrix_P() const
        {
            return m_perm.matrix_P<TT>();
        }
        
        Permute const &pivot_map() const
        {
            return m_perm;
        }
        
        // Get the packed LU matrix
        Mtx<TT> const matrix_LU() const { return m_matLU; }
        
        void deallocate() { m_matLU.deallocate(); m_perm.deallocate(); }
        
    private:
        Mtx<TT> m_matLU;
        Permute m_perm;
    };
    
}


#endif
