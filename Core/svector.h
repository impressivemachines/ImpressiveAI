//
//  svector.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/17/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_svector_h
#define Metaphor_svector_h

namespace im
{
    /*
    template <typename TT, int ROWS>
    class SVec
    {
    public:
        SVec() {}
        
        SVec(Vec<TT> const &v) { copy_from(v); }
        
        SVec(VecView<TT> const &vv) { copy_from(vv); }
        
        SVec(char const *pmatlabtext) { operator<<(pmatlabtext); }
        
        // Return the view
        VecView<TT> const &view() const { return VecView<TT>(ROWS, 1, m_data); }
        VecView<TT> view() { return VecView<TT>(ROWS, 1, m_data); }
        
        SVec &operator=(Vec<TT> const &v) { copy_from(v); return *this; }
        SVec &operator=(VecView<TT> const &vv) { copy_from(vv); return *this; }
        
        void copy_from(SVec const &sv) { memcpy(m_data, sv.m_data, ROWS * sizeof(TT)); }
        void copy_from(Vec<TT> const &v) { view().copy_from(v.view()); }
        void copy_from(VecView<TT> const &vv) { view().copy_from(vv); }
        
        // Copy over the source elements indicated by index list
        void copy_from(SVec const &src, std::vector<int> const &index_list)
        {
            view().copy_from(src.view(),index_list);
        }
        
        // Access information
        int rows() const { return ROWS; }
        int row_stride() const { return 1; }
        int count() const { return ROWS; }
        bool is_valid() const { return true; }
        
        // accessors
        TT const *ptr(int row) const
        {
            IM_DEBUG_ONLY_CHECK_BOUNDS(row, 0, ROWS);
            return m_data + row;
        }
        
        TT *ptr(int row)
        {
            IM_DEBUG_ONLY_CHECK_BOUNDS(row, 0, ROWS);
            return m_data + row;
        }
        
        TT const &operator()(int row) const { return *ptr(row); }
        TT &operator()(int row) { return *ptr(row); }
        
        TT const &at(int row) const { return *ptr(row); }
        TT &at(int row) { return *ptr(row); }
        
        TT const &index(int i) const { return *ptr(i); }
        TT &index(int i) { return *ptr(i); }
        
        // Paste over data from another view to a specific location.
        // This is done using block(...).copy_from(v.block(...)), but correctly handles mismatched
        // sizes and source/destination outside bounds, e.g. negative dstrowstart.
        
        void paste_from(SVec const &v, int dstrowstart)
        {
            view().paste_from(v.view(), dstrowstart);
        }
        
        void paste_from(Vec<TT> const &v, int dstrowstart)
        {
            view().paste_from(v.view(), dstrowstart);
        }
        
        // Fill with value
        void fill(TT const &f) { view().fill(f); }
        SVec &operator=(TT const &f) { view().fill(f); return *this; }
        
        // Assign the elements indicated by index list
        void fill(TT const &val, std::vector<int> const &index_list)
        {
            view().fill(val, index_list);
        }
        
        // Create from matlab format text field, e.g. "[ 1 2 3 4 5 6 ]" or "[3+6i; 4-5i]"
        SVec &operator<<(char const *ptext)
        {
            fill((TT)0);
            Mtx<TT> m(ptext);
            int size = std::max(m.rows(), m.cols());
            size = std::min(size, ROWS);
            if(m.cols()>m.rows())
                view().head(size).copy_from(m.view().row(0).head(size));
            else
                view().head(size).copy_from(m.view().col(0).head(size));
            return *this;
        }
        
        // Returns a square matrix with the vector elements as the diagonal
        Mtx<TT> diag_matrix() const
        {
            Mtx<TT> m(ROWS,ROWS);
            m = (TT)0;
            m.view().diag().copy_from(view());
            return m;
        }
        
        Mtx<TT> const matrix() const
        {
            return Mtx<TT>(ROWS,1,1,0,m_data);
        }
        
        Mtx<TT> matrix()
        {
            return Mtx<TT>(ROWS,1,1,0,m_data);
        }
        
        // Reshape the 1D matrix into a 2D matrix format
        
        // Row major
        Mtx<TT> const rm_2d_matrix(int colcount) const
        {
            IM_CHECK_ARGS(ROWS % colcount == 0);
            return Mtx<TT>(ROWS / colcount, colcount, colcount, 1, m_data);
        }
        
        Mtx<TT> rm_2d_matrix(int colcount)
        {
            IM_CHECK_ARGS(ROWS % colcount == 0);
            return Mtx<TT>(ROWS / colcount, colcount, colcount, 1, m_data);
        }
        
        // Column major
        Mtx<TT> const cm_2d_matrix(int rowcount) const
        {
            IM_CHECK_ARGS(ROWS % rowcount == 0);
            return Mtx<TT>(rowcount, ROWS / rowcount, 1, rowcount, m_data);
        }
        
        Mtx<TT> cm_2d_matrix(int rowcount)
        {
            IM_CHECK_ARGS(ROWS % rowcount == 0);
            return Mtx<TT>(rowcount, ROWS / rowcount, 1, rowcount, m_data);
        }
        
        // Copy over the source elements listed in the index list into consecutive rows
        void map(SVec const &src, std::vector<int> const &index_list)
        {
            view().map(src.view(), index_list);
        }
        
        
    private:
        TT m_data[ROWS];
    };
     */
}

#endif
