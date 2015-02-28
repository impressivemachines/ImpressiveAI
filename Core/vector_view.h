//
//  vector_view.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/8/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_vector_view_h
#define ImpressiveAI_vector_view_h

namespace im
{
    template <typename TT> class MtxView;
    
    template <typename TT>
    class VecView
    {
    public:
        typedef TT ValueType;
        
        VecView() : m_rows(0), m_rowstride(0), m_pdata(NULL) {}
        
        // Wrap vector-like memory
        VecView(int nrows, int rowstride, TT const *pdata)
        : m_rows(nrows), m_rowstride(rowstride), m_pdata(const_cast<TT *>(pdata))
        {
            IM_CHECK_LOWER_BOUNDS(nrows, 0);
        }

        // Wrap standard vector
        VecView(std::vector<TT> const &v)
        : m_rows((int)v.size()), m_rowstride(1), m_pdata(const_cast<TT *>(v.data()))
        {
        }
        
        void reset()
        {
            view_init(0, 0, NULL);
        }
        
        // Wrap vector-like memory
        // Returns a pointer to the first element after the vector in the source buffer
        TT const *wrap(int nrows, int rowstride, TT const *pdata)
        {
            view_init(nrows, rowstride, pdata);
            return pdata + nrows * rowstride;
        }
        
        // Wrap standard vector
        void wrap(const std::vector<TT> &v)
        {
            view_init((int)v.size(), 1, const_cast<TT *>(v.data()));
        }
        
        // Get info
        int rows() const { return m_rows; }
        int row_stride() const { return m_rowstride; }
        int count() const { return m_rows; }
        bool is_valid() const { return m_pdata!=NULL; }
        
        // Pointer to vector element
        TT const *ptr(int row) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_pdata);
            IM_DEBUG_ONLY_CHECK_BOUNDS(row, 0, m_rows);
            return m_pdata + row * m_rowstride;
        }
        
        TT *ptr(int row)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_pdata);
            IM_DEBUG_ONLY_CHECK_BOUNDS(row, 0, m_rows);
            return m_pdata + row * m_rowstride;
        }
        
        // Access the vector element
        TT const &operator()(int row) const { return *ptr(row); }
        TT &operator()(int row) { return *ptr(row); }
        
        TT const &at(int row) const { return *ptr(row); }
        TT &at(int row) { return *ptr(row); }
        
        // Get the raw data pointer
        // We do not check for null here
        TT const *ptr() const { return m_pdata; }
        TT *ptr() { return m_pdata; }
        
        // For consistency with matrix view
        TT const &index(int index) const { return at(index); }
        TT &index(int index) { return at(index); }
        
        // Copy over data from another view having the same size
        void copy_from(VecView const &vv)
        {
            IM_CHECK_VALID(vv);
            core_block_copy(*this, vv);
        }
        
        // Copy over the source elements indicated by index list
        void copy_from(VecView<TT> const &src, std::vector<int> const &index_list)
        {
            IM_CHECK_VALID(src);
            
            for(int i=0; i<index_list.size(); i++)
                at(index_list[i]) = src(index_list[i]);
        }
        
        // Paste over data from another view to a specific location.
        // This is done using block(...).copy_from(v.block(...)), but correctly handles mismatched
        // sizes and source/destination outside bounds, e.g. negative dstrowstart.
        void paste_from(VecView const &vv, int dstrowstart)
        {
            IM_CHECK_VALID(vv);
            
            int srcrowstart = 0;
            if(dstrowstart<0)
            {
                srcrowstart = -dstrowstart;
                dstrowstart = 0;
            }
            int rowcount = std::min(vv.rows() - srcrowstart, rows() - dstrowstart);
            if(rowcount>0)
                block(dstrowstart, rowcount).copy_from(vv.block(srcrowstart, rowcount));
        }
        
        // Set all elements to value
        void set(TT const &f)
        {
            core_block_fill(*this, f);
        }
        
        VecView &operator=(TT const &f)
        {
            set(f);
            return *this;
        }
        
        // Assign the elements indicated by index list
        void set(TT const &val, std::vector<int> const &index_list)
        {
            for(int i=0; i<index_list.size(); i++)
                at(index_list[i]) = val;
        }
        
        // Check equality
        bool operator==(VecView const &vv) const
        {
            IM_CHECK_VALID(vv);
            IM_CHECK_MATRIX_SIZES_MATCH(*this, vv);
            
            // Could be faster but probably only used for debugging
            int r;
            for(r=0; r<m_rows; r++)
                if(at(r)!=vv(r))
                    return false;
            return true;
        }
        
        bool operator!=(VecView const &vv) const
        {
            return !operator==(vv);
        }
        
        // Vector sub view
        VecView const block(int rowstart, int rowcount, int rowstep = 1) const
        {
            IM_CHECK_BOUNDS(rowstart, 0, m_rows);
            IM_CHECK_LOWER_BOUNDS(rowcount, 0);
            IM_CHECK_BOUNDS((rowstart + ( rowcount>0 ? (rowcount-1) * rowstep : 0)), 0, m_rows);
            
            return VecView(rowcount, rowstep * m_rowstride, ptr(rowstart));
        }
        
        VecView block(int rowstart, int rowcount, int rowstep = 1)
        {
            IM_CHECK_BOUNDS(rowstart, 0, m_rows);
            IM_CHECK_LOWER_BOUNDS(rowcount, 0);
            IM_CHECK_BOUNDS((rowstart + ( rowcount>0 ? (rowcount-1) * rowstep : 0)), 0, m_rows);
            
            return VecView(rowcount, rowstep * m_rowstride, ptr(rowstart));
        }
        
        // First few rows
        VecView const head(int nrows) const
        {
            nrows = std::min(nrows, m_rows);
            nrows = std::max(nrows, 0);
            
            return VecView(nrows, m_rowstride, ptr());
        }
        
        VecView head(int nrows)
        {
            nrows = std::min(nrows, m_rows);
            nrows = std::max(nrows, 0);

            return VecView(nrows, m_rowstride, ptr());
        }
        
        // Last few rows
        VecView const tail(int nrows) const
        {
            nrows = std::min(nrows, m_rows);
            nrows = std::max(nrows, 0);
            
            return VecView(nrows, m_rowstride, ptr(nrows > 0 ? m_rows - nrows : 0));
        }
        
        VecView tail(int nrows)
        {
            nrows = std::min(nrows, m_rows);
            nrows = std::max(nrows, 0);
            
            return VecView(nrows, m_rowstride, ptr(nrows > 0 ? m_rows - nrows : 0));
        }
        
        MtxView<TT> const matrix_view() const
        {
            return MtxView<TT>(m_rows, 1, m_rowstride, 0, ptr());
        }
        
        MtxView<TT> matrix_view()
        {
            return MtxView<TT>(m_rows, 1, m_rowstride, 0, ptr());
        }

        // Reshape the 1D matrix view into a 2D matrix format
        
        // Row major
        MtxView<TT> const rm_2d_matrix_view(int colcount) const
        {
            colcount = std::min(colcount, count());
            
            return MtxView<TT>(m_rows / colcount, colcount, m_rowstride * colcount, m_rowstride, ptr());
        }
        
        MtxView<TT> rm_2d_matrix_view(int colcount)
        {
            colcount = std::min(colcount, count());

            return MtxView<TT>(m_rows / colcount, colcount, m_rowstride * colcount, m_rowstride, ptr());
        }
        
        // Column major
        MtxView<TT> const cm_2d_matrix_view(int rowcount) const
        {
            rowcount = std::min(rowcount, count());
            
            return MtxView<TT>(rowcount, m_rows / rowcount, m_rowstride, m_rowstride * rowcount, ptr());
        }
        
        MtxView<TT> cm_2d_matrix_view(int rowcount)
        {
            rowcount = std::min(rowcount, count());
            
            return MtxView<TT>(rowcount, m_rows / rowcount, m_rowstride, m_rowstride * rowcount, ptr());
        }
        
        VecView const reverse() const
        {
            return VecView(m_rows, -m_rowstride, ptr(m_rows - 1));
        }
        
        VecView reverse()
        {
            return VecView(m_rows, -m_rowstride, ptr(m_rows - 1));
        }
        
        // Copy over the source elements listed in the index list into consecutive rows of *this.
        void map(VecView const &src, std::vector<int> const &index_list)
        {
            IM_CHECK_VALID(src);
            
            if(count()<1)
                return;
            
            int i;
            int row = 0;
            for(i=0; i<index_list.size(); i++)
            {
                at(row) = src(index_list[i]);
                row++;
                if(row==m_rows)
                    return;
            }
        }
        
        void print_size(FILE *fp = stdout) const
        {
            IM_CHECK_NULL(fp);
            
            fprintf(fp, "(%d)\n", m_rows);
        }
        
        // Print in Matlab format if known type
        void print(FILE *fp = stdout) const
        {
            IM_CHECK_NULL(fp);
    
            // Write on one line
            fprintf(fp, "[ ");
            for(int i=0; i<rows(); i++)
            {
                core_print_value(fp, at(i)); // prints "??" if the type is unknown
                if(i+1 < rows())
                    fprintf(fp, "; ");
            }
            fprintf(fp, " ]\n");

        }
        
        // Used to save as CSV etc.
        void write_text(FILE *fp, char delimiter = ' ') const
        {
            IM_CHECK_NULL(fp);
            
            for(int i=0; i<rows(); i++)
            {
                core_write_value(fp, at(i), delimiter);
                if(i+1 < rows())
                    fputc(delimiter, fp);
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
        
    protected:
        void view_init(int nrows, int rowstride, TT const *pdata)
        {
            IM_CHECK_LOWER_BOUNDS(nrows, 0);
            
            m_rows = nrows;
            m_rowstride = rowstride;
            m_pdata = const_cast<TT *>(pdata);
        }
        
    private:
        int m_rows;
        int m_rowstride;
        TT *m_pdata;
    };

    typedef VecView<uint8_t> VVb;
    typedef VecView<int16_t> VVs;
    typedef VecView<int> VVi;
    typedef VecView<float> VVf;
    typedef VecView<double> VVd;
    typedef VecView<im::Cf> VVcf;
    typedef VecView<im::Cd> VVcd;
}
#endif
