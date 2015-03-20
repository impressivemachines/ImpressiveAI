//
//  block_generic.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_generic_h
#define Metaphor_block_generic_h


namespace im
{
    template <typename TT> class MtxView;
    template <typename TT> class VecView;
    

    // Note that the source and destination matrix views can have very different orderings of memory
    // Optimization could be done for specific orderings and data types
    
    // core_block_fill - fill with value
    // core_block_copy - copy from one matrix to another
    // core_block_exchange - swap data between matrices
    // core_block_reshape - reshape by copying in row major order, copying as many elements as possible

    // Fill with value - works with any type
    template <typename TT> void core_block_fill(VecView<TT> vdst, TT const &val)
    {
        IM_CHECK_VALID(vdst);
        core_assign_loop(vdst.ptr(), vdst.row_stride(), val, vdst.rows());
    }
    
    template <typename TT> void core_block_fill(MtxView<TT> mdst, TT const &val)
    {
        IM_CHECK_VALID(mdst);
        
        int rows = mdst.rows();
        int cols = mdst.cols();
        
        if(rows>cols)
            for(int col = 0; col < cols; col++)
                core_assign_loop(mdst.ptr(0,col), mdst.row_stride(), val, rows);
        else
            for(int row = 0; row < rows; row++)
                core_assign_loop(mdst.ptr(row,0), mdst.col_stride(), val, cols);
    }

    // Copy over - works with any type
    
    template <typename TT> void core_block_copy(VecView<TT> vdst, VecView<TT> const &vsrc)
    {
        IM_CHECK_VALID(vdst);
        IM_CHECK_VALID(vsrc);
        IM_CHECK_VECTOR_SIZES_MATCH(vsrc, vdst);

        core_copy_loop(vdst.ptr(), vdst.row_stride(), vsrc.ptr(), vsrc.row_stride(), vdst.rows());
    }
    
    template <typename TT> void core_block_copy(MtxView<TT> mdst, MtxView<TT> const &msrc)
    {
        IM_CHECK_VALID(msrc);
        IM_CHECK_VALID(mdst);
        IM_CHECK_MATRIX_SIZES_MATCH(msrc,mdst);
        
        int rows = mdst.rows();
        int cols = mdst.cols();

        if(rows>cols)
            for(int col = 0; col < cols; col++)
                core_copy_loop(mdst.ptr(0,col), mdst.row_stride(), msrc.ptr(0,col), msrc.row_stride(), rows);
        else
            for(int row = 0; row < rows; row++)
                core_copy_loop(mdst.ptr(row,0), mdst.col_stride(), msrc.ptr(row,0), msrc.col_stride(), cols);
    }

    // Exchange values - works with any type
    
    template <typename TT> void core_block_exchange(VecView<TT> v1, VecView<TT> v2)
    {
        IM_CHECK_VALID(v1);
        IM_CHECK_VALID(v2);
        IM_CHECK_VECTOR_SIZES_MATCH(v1, v2);
        
        core_exchange_loop(v1.ptr(), v1.row_stride(), v2.ptr(), v2.row_stride(), v1.rows());
    }
    
    template <typename TT> void core_block_exchange(MtxView<TT> m1, MtxView<TT> m2)
    {
        IM_CHECK_VALID(m1);
        IM_CHECK_VALID(m2);
        IM_CHECK_MATRIX_SIZES_MATCH(m1,m2);
        
        int rows = m1.rows();
        int cols = m1.cols();
        
        if(rows>cols)
            for(int col = 0; col < cols; col++)
                core_exchange_loop(m1.ptr(0,col), m1.row_stride(), m2.ptr(0,col), m2.row_stride(), rows);
        else
            for(int row = 0; row < rows; row++)
                core_exchange_loop(m1.ptr(row,0), m1.col_stride(), m2.ptr(row,0), m2.col_stride(), cols);
    }

    // Reshape by copying in row major raster order - works with any type
    template <typename TT> void core_block_reshape(MtxView<TT> mdst, MtxView<TT> const &msrc)
    {
        IM_CHECK_VALID(msrc);
        IM_CHECK_VALID(mdst);
        
        int srcrows = msrc.rows();
        int srccols = msrc.cols();
        int dstrows = mdst.rows();
        int dstcols = mdst.cols();
        
        if(msrc.count() < mdst.count())
        {
            int n = msrc.count() / dstcols;
            mdst.block(n,0,dstrows-n,dstcols) = (TT)0;
        }
        
        int dststride = mdst.col_stride();
        
        int sr = 0;
        int sc = 0;
        
        for(int row = 0; row<dstrows; row++)
        {
            TT *pdst = mdst.ptr(row,0);
            for(int col = 0; col<dstcols; col++, pdst += dststride)
            {
                pdst[0] = msrc(sr,sc);
                
                sc++;
                if(sc<srccols)
                    continue;
                
                sc = 0;
                sr++;
                
                if(sr<srcrows)
                    continue;
                
                return;
            }
        }
    }
    
}



#endif
