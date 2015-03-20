//
//  decomp_qr.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_decomp_qr_h
#define Metaphor_decomp_qr_h

namespace im
{
    // Factorise a general M x N matrix A into
    // A P = Q R
    // Where Q is orthogonal M x M and R is upper (right) triangular M x N
    // P is the pivot matrix - if use_pivoting is false then this is set to the identity.
    // If A is rank deficient with rank r<n, then P is used to ensure that the lower n-r
    // rows of R are zero and the first r columns of Q are an orthonormal basis for A.
    //
    // Alternatively, the thin version of Q and R can be computed. In this case,
    // Q is M x k and R is k x N, where k is min(M,N). This is better for M>>N. The k
    // columns of Q are orthornomal and for M>=N, R is square and upper right triangular.
    
    template <typename TT>
    class MatrixDecompQR
    {
    public:
        MatrixDecompQR() {}
        MatrixDecompQR(MtxView<TT> const &mavA) { compute(mavA); }
        
        // Pivoting is useful when A is nearly rank deficient
        void compute(MtxView<TT> const &mavA, bool use_pivoting = true);
        
        // Solve AX = Y using Q R P^T factorization (for square A only)
        Mtx<TT> solve(MtxView<TT> const &mavy) const;
        Vec<TT> solve(VecView<TT> const &vvy) const;
        
        // Solves the overdetermined system AX=Y for A.rows() > A.cols() using least squares
        Mtx<TT> solve_least_squares(MtxView<TT> const &mavy) const;
        
        // Compute the inverse
        Mtx<TT> inverse() const;
        
        // Returns the determinant
        TT det() const
        {
            if(m_matQR.rows()!=m_matQR.cols())
                return (TT)0;
            
            TT det = (TT)m_perm.sign();
            for(int i=0; i<m_matQR.rows(); i++)
                det *= m_matQR(i,i);
            
            return det;
        }
        
        // Returns the log of the abs(determinant)
        TT log_det() const
        {
            if(m_matQR.rows()!=m_matQR.cols())
                IM_THROW_NO_SOLUTION;
            
            TT logdet = (TT)0;
            for(int i=0; i<m_matQR.rows(); i++)
                logdet += std::log(std::abs(m_matQR(i,i)));
            
            return logdet;
        }
        
        // Returns the sign of the determinant, or 0 if singular
        TT sign_det() const
        {
            if(m_matQR.rows()!=m_matQR.cols())
                return (TT)0;
            
            TT s = (TT)m_perm.sign();
            for(int i=0; i<m_matQR.rows(); i++)
            {
                TT val = m_matQR(i,i);
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
        
        Mtx<TT> const matrix_Q()
        {
            unpack_QR();
            return m_matQ;
        }
        
        Mtx<TT> const matrix_R()
        {
            unpack_QR();
            return m_matR;
        }
        
        Mtx<TT> const matrix_Q_thin()
        {
            unpack_QR_thin();
            return m_matQ;
        }
        
        Mtx<TT> const matrix_R_thin()
        {
            unpack_QR_thin();
            return m_matR;
        }
        
        // permutation matrix - this will be the identity if use_pivoting is false
        Mtx<TT> matrix_P() const
        {
            return m_perm.matrix_P<TT>().t();
        }
        
        // Compute Q^T X efficiently
        Mtx<TT> matrix_QTX(MtxView<TT> const &mavX) const
        {
            Mtx<TT> matY(mavX.rows(), mavX.cols());
            matY.copy_from(mavX);
            matrix_QTX_in_place(matY.view());
            return matY;
        }
        
         // Compute Q^T X efficiently
        void matrix_QTX_in_place(MtxView<TT> mavX) const;

        void deallocate()
        {
            m_matQR.deallocate();
            m_matQ.deallocate();
            m_matR.deallocate();
            m_vTau.deallocate();
            m_perm.deallocate();
        }
        
    private:
        void unpack_QR();
        void unpack_QR_thin();
        
        Mtx<TT> m_matQR; // Q stored as householder transformations in lower tri part, R in upper tri + diag
        Mtx<TT> m_matQ, m_matR; // Stores unpacked version if required
        bool m_unpacked;
        bool m_unpacked_thin;
        Vec<TT> m_vTau;
        Permute m_perm;
    };
}

#endif
