//
//  convolve.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/3/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_convolve_h
#define Metaphor_convolve_h

namespace im
{
    // This class gets a set of standard size rectangles to process from a source rectangular
    // region. The next() call gets the next sub rectangle to process or returns false if there are
    // none left. The block rect returned is located with the bounds of the source rect, and is
    // normally of size block_size, but will be smaller for blocks overlapping the left or bottom edge.
    
    class RectSubdivider
    {
    public:
        RectSubdivider() {}
        
        // Source is the (larger) containing rectangle, and block_size is the maximum size of a block
        // to return from next().
        void init(MtxRect const &source, MtxSize const &block_size)
        {
            m_source_rect = source;
            m_block_size = block_size;
            reset();
        }
        
        void reset() { m_loc = MtxLoc(0,0); }
        
        bool is_active() const { return m_loc.row < m_source_rect.size.rows; }
        
        // Get the next block or return false if there are none left to process.
        bool next(MtxRect &block)
        {
            if(!is_active())
                return false;
            
            block.origin.row = m_source_rect.origin.row + m_loc.row;
            block.origin.col = m_source_rect.origin.col + m_loc.col;
            block.size.rows = std::min(m_block_size.rows, m_source_rect.size.rows - m_loc.row);
            block.size.cols = std::min(m_block_size.cols, m_source_rect.size.cols - m_loc.col);
            
            m_loc.col += m_block_size.cols;
            if(m_loc.col >= m_source_rect.size.cols)
            {
                m_loc.col = 0;
                m_loc.row += m_block_size.rows;
            }

            return true;
        }
        
    private:
        MtxRect m_source_rect; // input rect
        MtxSize m_block_size; // size of block
        MtxLoc m_loc; // current top left location of block relative to rect origin
    };
    
    // This class gets typically fixed sized rectangular blocks from a given view.
    // If the rect of the desired block is partly or fully outside the source matrix, then
    // the class does the correct padding by copying into an internal block store and
    // returning a reference to that instead.
    
    template <typename TT>
    class MatrixPadder
    {
    public:
        MatrixPadder() {}
        MatrixPadder(MtxView<TT> const &mv) : m_source(mv) {}
        
        void init(MtxView<TT> const &mv, int block_size = 0)
        {
            m_source = mv;
            m_block_store.resize(block_size, block_size);
        }
        
        MtxView<TT> get(MtxRect const &rct, PadMode pad);
        
        void deallocate() { m_block_store.deallocate(); }
        
    private:
        void pad_zero(MtxRect const &rct);
        void pad_extend(MtxRect const &rct);
        void pad_wrap(MtxRect const &rct);
        void pad_reflect(MtxRect const &rct);
        
        MtxView<TT> m_source;
        Mtx<TT> m_block_store;
    };

    // General 2D convolution
    // The kernel center indicates the zero origin of the kernel matrix
    // The start location indicates the position in the input matrix corresponding to the top left of the output matrix
    // output_pre_weight is the multiplier for the existing output before adding the convolution results, e.g. 0 for
    // normal convolution, and 1 for adding the result to the output.
    // drows and dcols are the output decimation factors.
    // Tune blocksize for best performace.
    template <typename TT>
    void core_convolve(MtxView<TT> mavout, MtxView<TT> const &mavkernel, MtxView<TT> const &mavin,
                              PadMode pad, MtxLoc kernel_center, MtxLoc start_loc, TT output_pre_weight, int drows, int dcols, int blocksize = 32);
    
    // Specialization for 1D
    // The kernel center indicates the zero origin of the kernel matrix
    // The start location indicates the position in the input matrix corresponding to the first sample of the output matrix
    // output_pre_weight is the multiplier for the existing output before adding the convolution results, e.g. 0 for
    // normal convolution, and 1 for adding the result to the output.
    template <typename TT>
    void core_convolve(VecView<TT> vvout, VecView<TT> const &vvkernel, VecView<TT> const &vvin,
                              PadMode pad, int kernel_center, int start_loc, TT output_pre_weight, int drows);
    
}

#endif
