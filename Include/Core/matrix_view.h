//
//  matrix_view.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 6/17/14.
//  Copyright (c) 2014 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpAILibrary_matrix_view_h
#define ImpAILibrary_matrix_view_h

namespace im
{
    // The matrix view is a light weight generic way for referencing existing 1D or 2D memory in
    // matrix-like format.
    // It has a number of sub-view and re-ordering
    // methods, none of which do any data copying. Note that the view does not manage the underlying memory.
    // Use Mtx if you want to have a resizable object that behaves like MtxView and manages its own memory.
    // We decided to make one object which works with both const and non-const buffers.
    // If you want the data to be const, then construct a MtxView const.
    // Note that operations do not check if the same memory is used for source and destination.
    // MtxView works on all types TT
    // Makes use of BlockCopy, BlockReshape from block_generic.h
    
    template <typename TT> class MtxViewIterator;
    
    template <typename TT>
    class MtxView
    {
    public:
        typedef TT ValueType;
        
        MtxView() : m_rows(0), m_cols(0), m_rowstride(0), m_colstride(0), m_pdata(NULL) {}
        
        // Wrap matrix-like memory
        MtxView(int nrows, int ncols, int rowstride, int colstride, TT const *pdata)
        : m_rows(nrows), m_cols(ncols), m_rowstride(rowstride), m_colstride(colstride), m_pdata(const_cast<TT *>(pdata))
        {
            IM_CHECK_LOWER_BOUNDS(nrows, 0);
            IM_CHECK_LOWER_BOUNDS(ncols, 0);
        }
        
        // Wrap standard vector
        MtxView(std::vector<TT> const &v)
        : m_rows((int)v.size()), m_cols(1), m_rowstride(1), m_colstride(0), m_pdata(const_cast<TT *>(v.data()))
        {
        }
        
        void reset()
        {
            view_init(0, 0, 0, 0, NULL);
        }
        
        // Wrap matrix-like memory
        // Returns a pointer to the first element after the matrix in the source buffer
        TT const *wrap(int nrows, int ncols, int rowstride, int colstride, TT const *pdata)
        {
            view_init(nrows, ncols, rowstride, colstride, pdata);
            return pdata + nrows * rowstride;
        }
        
        // Wrap standard vector
        void wrap(const std::vector<TT> &v)
        {
            view_init((int)v.size(), 1, 1, 0, const_cast<TT *>(v.data()));
        }
        
        // Get info
        int rows() const { return m_rows; }
        int cols() const { return m_cols; }
        int row_stride() const { return m_rowstride; }
        int col_stride() const { return m_colstride; }
        int count() const { return m_rows * m_cols; }
        bool is_valid() const { return m_pdata!=NULL; }
        
        // Pointer to matrix element
        TT const *ptr(int row, int col) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_pdata);
            IM_DEBUG_ONLY_CHECK_BOUNDS(row, 0, m_rows);
            IM_DEBUG_ONLY_CHECK_BOUNDS(col, 0, m_cols);
            return m_pdata + row * m_rowstride + col * m_colstride;
        }
        
        TT *ptr(int row, int col)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_pdata);
            IM_DEBUG_ONLY_CHECK_BOUNDS(row, 0, m_rows);
            IM_DEBUG_ONLY_CHECK_BOUNDS(col, 0, m_cols);
            return m_pdata + row * m_rowstride + col * m_colstride;
        }
        
        // Access the matrix element
        TT const &operator()(int row, int col) const { return *ptr(row, col); }
        TT &operator()(int row, int col) { return *ptr(row, col); }
        
        TT const &at(int row, int col) const { return *ptr(row, col); }
        TT &at(int row, int col) { return *ptr(row, col); }
        
        // Get the raw data pointer
        // We do not check for null here
        TT const *ptr() const { return m_pdata; }
        TT *ptr() { return m_pdata; }
        
        // Index as linear array regardless of shape
        TT const &index(int i) const { return at(i / m_cols, i % m_cols); }
        TT &index(int i) { return at(i / m_cols, i % m_cols); }
        
        // Copy over data from another view having the same size
        void copy_from(MtxView const &mav)
        {
            IM_CHECK_VALID(mav);
            
            core_block_copy(*this, mav);
        }
        
        // Copy over the source elements indicated by index list
        void copy_from(MtxView<TT> const &src, std::vector<MtxLoc> const &index_list)
        {
            IM_CHECK_VALID(src);
            
            for(int i=0; i<index_list.size(); i++)
            {
                const MtxLoc &loc = index_list[i];
                at(loc.row, loc.col) = src(loc.row, loc.col);
            }
        }
        
        // Copy data in raster order without regard to the relative shapes.
        // Any destination elements which have no source data are set to zero.
        void reshape_from(MtxView const &mav)
        {
            IM_CHECK_VALID(mav);
            
            if(rows()>0 && cols()>0)
                core_block_reshape(*this, mav);
        }
        
        // Paste over data from another view to a specific location.
        // This is done using Block(...).CopyFrom(mv.Block(...)), but correctly handles mismatched
        // sizes and source/destination outside bounds, e.g. negative dstrowstart.
        void paste_from(MtxView const &mav, int dstrowstart, int dstcolstart)
        {
            IM_CHECK_VALID(mav);
            
            int srcrowstart = 0;
            if(dstrowstart<0)
            {
                srcrowstart = -dstrowstart;
                dstrowstart = 0;
            }
            int rowcount = std::min(mav.rows() - srcrowstart, rows() - dstrowstart);
            int srccolstart = 0;
            if(dstcolstart<0)
            {
                srccolstart = -dstcolstart;
                dstcolstart = 0;
            }
            int colcount = std::min(mav.cols() - srccolstart, cols() - dstcolstart);
            if(rowcount>0 && colcount>0)
            {
                block(dstrowstart, dstcolstart, rowcount, colcount).copy_from(mav.block(srcrowstart, srccolstart, rowcount, colcount));
            }
        }
        
        // Set all elements to value
        void fill(TT const &f)
        {
            core_block_fill(*this, f);
        }
        
        MtxView &operator=(TT const &f)
        {
            fill(f);
            return *this;
        }
        
        
        // Assign the elements indicated by index list
        void fill(TT const &val, std::vector<MtxLoc> const &index_list)
        {
            for(int i=0; i<index_list.size(); i++)
            {
                const MtxLoc &loc = index_list[i];
                at(loc.row, loc.col) = val;
            }
        }
        
        // Check equality
        bool operator==(MtxView const &mav) const
        {
            IM_CHECK_VALID(mav);
            IM_CHECK_MATRIX_SIZES_MATCH(*this, mav);
            
            // Could be faster but probably only used for debugging
            int r,c;
            for(r=0; r<m_rows; r++)
                for(c=0; c<m_cols; c++)
                    if(at(r,c)!=mav(r,c))
                        return false;
            return true;
        }
        
        bool operator!=(MtxView const &mav) const
        {
            return !operator==(mav);
        }
        
        // Matrix sub view with optional decimation and/or reflection
        MtxView const block(int rowstart, int colstart,
                                   int rowcount, int colcount,
                                   int rowstep = 1, int colstep = 1) const
        {
            IM_CHECK_BOUNDS(rowstart, 0, m_rows);
            IM_CHECK_BOUNDS(colstart, 0, m_cols);
            IM_CHECK_LOWER_BOUNDS(rowcount, 0);
            IM_CHECK_LOWER_BOUNDS(colcount, 0);
            IM_CHECK_BOUNDS((rowstart + ( rowcount>0 ? (rowcount-1) * rowstep : 0)), 0, m_rows);
            IM_CHECK_BOUNDS((colstart + ( colcount>0 ? (colcount-1) * colstep : 0)), 0, m_cols);
            
            return MtxView(rowcount, colcount, rowstep * m_rowstride, colstep * m_colstride, ptr(rowstart, colstart));
        }
        
        MtxView block(int rowstart, int colstart,
                            int rowcount, int colcount,
                            int rowstep = 1, int colstep = 1)
        {
            IM_CHECK_BOUNDS(rowstart, 0, m_rows);
            IM_CHECK_BOUNDS(colstart, 0, m_cols);
            IM_CHECK_LOWER_BOUNDS(rowcount, 0);
            IM_CHECK_LOWER_BOUNDS(colcount, 0);
            IM_CHECK_BOUNDS((rowstart + ( rowcount>0 ? (rowcount-1) * rowstep : 0)), 0, m_rows);
            IM_CHECK_BOUNDS((colstart + ( colcount>0 ? (colcount-1) * colstep : 0)), 0, m_cols);
            
            return MtxView(rowcount, colcount, rowstep * m_rowstride, colstep * m_colstride, ptr(rowstart, colstart));
        }
        
        MtxView const block(MtxRect const &rct) const { return block(rct.origin.row, rct.origin.col, rct.size.rows, rct.size.cols); }
        MtxView block(MtxRect const &rct) { return block(rct.origin.row, rct.origin.col, rct.size.rows, rct.size.cols); }
        
        // Return a row as a vector
        VecView<TT> const row(int r) const
        {
            IM_CHECK_BOUNDS(r, 0, m_rows);
       
            return VecView<TT>(m_cols, m_colstride, ptr(r, 0));
        }
        
        VecView<TT> row(int r)
        {
            IM_CHECK_BOUNDS(r, 0, m_rows);
            
            return VecView<TT>(m_cols, m_colstride, ptr(r, 0));
        }
        
        // Return a column as a vector
        VecView<TT> const col(int c) const
        {
            IM_CHECK_BOUNDS(c, 0, m_cols);
            return VecView<TT>(m_rows, m_rowstride, ptr(0, c));
        }
        
        VecView<TT> col(int c)
        {
            IM_CHECK_BOUNDS(c, 0, m_cols);
            return VecView<TT>(m_rows, m_rowstride, ptr(0, c));
        }
        
        // Return the diagonal as a vector
        // Positive offset gets the diagonal starting from (0,offset)
        // Negative offset gets the diagonal starting from (-offset,0)
        VecView<TT> const diag(int offset = 0) const
        {
            if(offset>=0)
            {
                int size = std::min(rows(), cols()-offset);
                size = std::max(size,0);
                return VecView<TT>(size, row_stride() + col_stride(), ptr(0, offset));
            }
            else
            {
                int size = std::min(rows()+offset, cols());
                size = std::max(size,0);
                return VecView<TT>(size, row_stride() + col_stride(), ptr(-offset, 0));
            }
        }
        
        VecView<TT> diag(int offset = 0)
        {
            if(offset>=0)
            {
                int size = std::min(rows(), cols()-offset);
                size = std::max(size,0);
                return VecView<TT>(size, row_stride() + col_stride(), ptr(0, offset));
            }
            else
            {
                int size = std::min(rows()+offset, cols());
                size = std::max(size,0);
                return VecView<TT>(size, row_stride() + col_stride(), ptr(-offset, 0));
            }
        }
        
        // Return the transpose matrix or vector
        MtxView const t() const
        {
            return MtxView(m_cols, m_rows, m_colstride, m_rowstride, ptr());
        }
        
        MtxView t()
        {
            return MtxView(m_cols, m_rows, m_colstride, m_rowstride, ptr());
        }
        
        // Return this matrix flipped horizontally
        MtxView const reverse_cols() const
        {
            return MtxView(m_rows, m_cols, m_rowstride, -m_colstride, ptr(0, m_cols - 1));
        }
        
        MtxView reverse_cols()
        {
            return MtxView(m_rows, m_cols, m_rowstride, -m_colstride, ptr(0, m_cols - 1));
        }
        
        // Return this matrix flipped vertically
        MtxView const reverse_rows() const
        {
            return MtxView(m_rows, m_cols, -m_rowstride, m_colstride, ptr(m_rows - 1, 0));
        }
        
        MtxView reverse_rows()
        {
            return MtxView(m_rows, m_cols, -m_rowstride, m_colstride, ptr(m_rows - 1, 0));
        }
        
        MtxView const rotate_90cw() const
        {
            return t().reverse_cols();
        }
        
        MtxView rotate_90cw()
        {
            return t().reverse_cols();
        }
        
        MtxView const rotate_90ccw() const
        {
            return t().reverse_rows();
        }
        
        MtxView rotate_90ccw()
        {
            return t().reverse_rows();
        }
        
        MtxView const rotate_180() const
        {
            return reverse_cols().reverse_rows();
        }
        
        MtxView rotate_180()
        {
            return reverse_cols().reverse_rows();
        }

        // Copy over the source elements listed in the index list into consecutive cols/rows of *this.
        void map(MtxView const &src, std::vector<MtxLoc> const &index_list)
        {
            IM_CHECK_VALID(src);
            
            if(count()<1)
                return;
            
            int i;
            int row = 0;
            int col = 0;
            for(i=0; i<index_list.size(); i++)
            {
                const MtxLoc &loc = index_list[i];
                at(row,col) = src(loc.row, loc.col);
                col++;
                if(col<m_cols)
                    continue;
                col = 0;
                row++;
                if(row<m_rows)
                    continue;
                return;
            }
        }
        
        // Copy over the source rows listed in the index list into consecutive rows of *this.
        void map_rows(MtxView const &src, std::vector<int> const &index_list)
        {
            IM_CHECK_VALID(src);
            IM_CHECK_ARGS(src.cols()==cols());
            IM_CHECK_ARGS(index_list.size()<=rows());
            
            for(int r=0; r<index_list.size(); r++)
                core_block_copy(row(r), src.row(index_list[r]));
        }
       
        // Copy over the source cols listed in the index list into consecutive cols of *this.
        void map_cols(MtxView const &src, std::vector<int> const &index_list)
        {
            IM_CHECK_VALID(src);
            IM_CHECK_ARGS(src.rows()==rows());
            IM_CHECK_ARGS(index_list.size()<=cols());
            
            for(int c=0; c<index_list.size(); c++)
                core_block_copy(col(c), src.col(index_list[c]));
        }
        
        void print_size(bool cr = true, FILE *fp = stdout) const
        {
            IM_CHECK_NULL(fp);
            
            fprintf(fp, "(%d x %d)", m_rows, m_cols);
            if(cr)
                fputc('\n', fp);
        }
        
        // Print in Matlab format if known type
        void print(FILE *fp = stdout) const
        {
            IM_CHECK_NULL(fp);
            
            int nrows = rows();
            int ncols = cols();
            
            if(ncols==1)
            {
                // Write on one line
                fprintf(fp, "[ ");
                for(int i=0; i<nrows; i++)
                {
                    core_print_value(fp, at(i,0)); // prints "??" if the type is unknown
                    if(i+1 < nrows)
                        fprintf(fp, "; ");
                }
                fprintf(fp, " ]\n");
            }
            else if(nrows==1)
            {
                // Write on one line
                fprintf(fp, "[ ");
                for(int i=0; i<ncols; i++)
                {
                    core_print_value(fp, at(0,i));
                    if(i+1 < ncols)
                        fprintf(fp, "  ");
                }
                fprintf(fp, " ]\n");
            }
            else
            {
                fprintf(fp, "[");
                for(int r=0; r<nrows; r++)
                {
                    if(r!=0)
                        fprintf(fp, " ");
                    for(int c=0; c<ncols; c++)
                    {
                        fprintf(fp, " ");
                        core_print_value(fp, at(r,c));
                    }
                    if(r+1 < nrows)
                        fprintf(fp, ";\n");
                    else
                        fprintf(fp, " ]\n");
                }
            }
        }

        // Used to save as CSV etc.
        void write_text(FILE *fp, char delimiter = ' ') const
        {
            IM_CHECK_NULL(fp);
            
            for(int row=0; row<rows(); row++)
            {
                for(int col=0; col<cols(); col++)
                {
                    core_write_value(fp, at(row,col), delimiter);
                    if(col<cols()-1)
                        fputc(delimiter, fp);
                }
                
                fputc('\n', fp);
            }
            
            if(ferror(fp))
                IM_THROW_FILE_ERROR;
        }
        
        // Save as a file using WriteText
        void save_text(char const *pfile, char delimiter = ' ') const
        {
            IM_CHECK_NULL(pfile);
            
            FILE *fp = fopen(pfile, "w");
            if(fp==NULL)
                IM_THROW_FILE_ERROR;
            
            write_text(fp, delimiter);
            
            fclose(fp);
        }
        
        MtxViewIterator<TT> iterator_begin()
        {
            MtxViewIterator<TT> it(*this);
            return it;
        }
        
        MtxViewIterator<TT> iterator_end()
        {
            MtxViewIterator<TT> it(*this);
            it.goto_end();
            return it;
        }
        
    protected:
        void view_init(int nrows, int ncols, int rowstride, int colstride, TT const *pdata)
        {
            IM_CHECK_LOWER_BOUNDS(nrows, 0);
            IM_CHECK_LOWER_BOUNDS(ncols, 0);
            
            m_rows = nrows;
            m_cols = ncols;
            m_rowstride = rowstride;
            m_colstride = colstride;
            m_pdata = const_cast<TT *>(pdata);
        }
        
    private:
        int m_rows, m_cols;
        int m_rowstride, m_colstride;
        TT *m_pdata;
    };
    
    template <typename TT>
    class MtxViewIterator : public std::iterator<std::forward_iterator_tag, TT>
    {
    public:
        MtxViewIterator() : row(0), col(0) {}
        MtxViewIterator(MtxView<TT> const &mav) : row(0), col(0), m_mav(mav) {}
        
        TT *ptr() { return m_mav.ptr(row,col); }
        
        TT *ptr(int drow, int dcol)
        {
            int r = row + drow;
            int c = col + dcol;
            
            // Clip to bounds of matrix size
            r = std::max(r, 0);
            r = std::min(r, m_mav.rows()-1);
            c = std::max(c, 0);
            c = std::min(c, m_mav.cols()-1);
            
            return m_mav.ptr(r, c);
        }
        
        TT &at() { return *ptr(); }
        
        TT &at(int drow, int dcol) { return *ptr(drow, dcol); }
        
        TT &operator*() { return at(); }
        
        void goto_start()
        {
            row = 0;
            col = 0;
        }
        
        void goto_end()
        {
            row = m_mav.rows();
            col = 0;
        }

        bool is_active() const { return row<m_mav.rows(); }
        
        bool next()
        {
            if(!is_active())
                return false;
            
            if(col==m_mav.cols()-1)
            {
                row++;
                col = 0;
                if(!is_active())
                    return false;
            }
            else
                col++;
            
            return true;
        }
        
        MtxView<TT> view() const { return m_mav; }
        
        MtxViewIterator operator++(int dummy) { MtxViewIterator<TT> tmp(*this); next(); return tmp; } // postfix increment iter++
        MtxViewIterator &operator++() { next(); return *this; } // prefix increment ++iter
        
        bool operator==(MtxViewIterator const &other) { return row==other.row && col==other.col; }
        bool operator!=(MtxViewIterator const &other) { return !operator==(other); }
        
    public:
        int row, col; // you can simply assign these to go to a specific location
        
    private:
        MtxView<TT> m_mav;
    };
    
    typedef MtxView<uint8_t> MVb;
    typedef MtxView<int16_t> MVs;
    typedef MtxView<int> MVi;
    typedef MtxView<float> MVf;
    typedef MtxView<double> MVd;
    typedef MtxView<im::Cf> MVcf;
    typedef MtxView<im::Cd> MVcd;
    
    typedef MtxViewIterator<uint8_t> MVIb;
    typedef MtxViewIterator<int16_t> MVIs;
    typedef MtxViewIterator<int> MVIi;
    typedef MtxViewIterator<float> MVIf;
    typedef MtxViewIterator<double> MVId;
    typedef MtxViewIterator<im::Cf> MVIcf;
    typedef MtxViewIterator<im::Cd> MVIcd;
}



#endif
