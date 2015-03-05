//
//  permute.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_permute_h
#define Metaphor_permute_h

namespace im
{
    // The permute class maintains a permutation matrix, based on swaps of indices in a vector
    // This can then be used to multiply any given matrix in order to permute rows or columns
    
    // The P matrix is a square matrix that starts out as the identity
    // When you call swap, it swaps the indicated rows of P
    
    class Permute
    {
    public:
        Permute() {}
        Permute(int size) { resize(size); }
        
        void resize(int size) { m_permute_map.resize(size); init(); }
        
        // Set back to identity
        void init()
        {
            for(int i=0; i<m_permute_map.size(); i++)
                m_permute_map[i] = i;
            m_sign = 1;
        }
        
        // swap rows of P, cols of P^T
        void swap(int i, int j)
        {
            if(i==j)
                return;
            
            IM_DEBUG_ONLY_CHECK_BOUNDS(i, 0, m_permute_map.size());
            IM_DEBUG_ONLY_CHECK_BOUNDS(j, 0, m_permute_map.size());
            
            std::swap(m_permute_map[i], m_permute_map[j]);
            m_sign = -m_sign;
        }
        
        int const &operator()(int i) const
        {
            IM_DEBUG_ONLY_CHECK_BOUNDS(i, 0, m_permute_map.size());
            return m_permute_map[i];
        }
        
        int &operator()(int i)
        {
            IM_DEBUG_ONLY_CHECK_BOUNDS(i, 0, m_permute_map.size());
            return m_permute_map[i];
        }
        
        // positive for even numbers of swaps
        int sign() const { return m_sign; }
        
        // Get the permute map as a matrix
        template <typename TT>
        Mtx<TT> matrix_P() const
        {
            int size = (int)m_permute_map.size();
            Mtx<TT> mat(size, size);
            mat = (TT)0;
            for(int i=0; i<size; i++)
                mat(i, m_permute_map[i]) = (TT)1;
            return mat;
        }
        
        // Permutes the rows of X: X = PX
        template <typename TT>
        void matrix_PX_in_place(MtxView<TT> mavX) const;
        
        // Un-permutes the rows of X: X = P^T X
        template <typename TT>
        void matrix_PTX_in_place(MtxView<TT> mavX) const;
        
        void deallocate() { m_permute_map.clear(); }
        
    private:
        std::vector<int> m_permute_map;
        int m_sign;
    };
}


#endif
