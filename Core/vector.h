//
//  vector.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/19/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_vector_h
#define ImpressiveAI_vector_h

namespace im
{
    template <typename TT> class Mtx;
    
    // While VecView is just a view on memory, Vec also manages memory.
    // Vec also defines a number of additional operations such as vector operators +* etc.
    // Vec is only implemented for float, double, complex float and complex double.
    // Vec can wrap external memory by initializing it from a VecView.
    
    template <typename TT>
    class Vec
    {
    private:
        VecView<TT> m_view;
        std::shared_ptr<std::vector<TT>> m_mem;
        
    protected:
        friend class Mtx<TT>;
        
        Vec(std::shared_ptr<std::vector<TT>> const &mem, VecView<TT> const &vv) : m_mem(mem), m_view(vv) {}
        
    public:
        typedef TT ValueType;
        
        Vec() { resize(0); }
        
        // Construct a vector
        Vec(int nrows) { resize(nrows); }
        
        // Construct using << operator from matlab formatted text string representation, eg. "[ 1 2 3 4 ]"
        Vec(char const *pmatlabtext) { *this << pmatlabtext; }
        
        // Wrap the given view - note that this object then manages no memory and that there is no data copy
        Vec(VecView<TT> const &vv) { m_mem.reset(); m_view = vv; }
        
        Vec(int rows, int rowstride, TT const *pdata)
        {
            m_mem.reset();
            m_view.wrap(rows, rowstride, pdata);
        }
        
        TT const *wrap(int nrows, int rowstride, TT const *pdata)
        {
            m_mem.reset();
            return m_view.wrap(nrows, rowstride, pdata);
        }
        
        void wrap(const std::vector<TT> &v)
        {
            m_mem.reset();
            m_view.wrap(v);
        }
        
        // Resizes the matrix and loses any data
        void resize(int nrows);
        
        // Wrap the given view - note that this object then manages no memory and that there is no data copy
        Vec &operator=(VecView<TT> const &vv) { m_mem.reset(); m_view = vv; return *this; }
        
        // Return the view
        VecView<TT> const &view() const { return m_view; }
        VecView<TT> view() { return m_view; }
        
        void deallocate() { m_view.reset(); m_mem.reset(); }
        
        // Stops sharing the vector by creating own private copy and guarentees packed layout
        void stop_sharing();
        
        // Ensures vector has the requested capacity non-destructively
        // This allows future resize operations to proceed without allocation overhead
        void reserve(int nrows)
        {
            IM_CHECK_LOWER_BOUNDS(nrows, 0);
            
            if(m_mem)
                m_mem->reserve(nrows);
        }
        
        // Add an new row at the bottom of the vector and sets its element to the given value
        void push_back(TT const &f);
        
        // Remove a row from the bottom of the matrix
        void pop_back();
        
        // Access information
        int rows() const { return m_view.rows(); }
        int row_stride() const { return m_view.row_stride(); }
        int count() const { return m_view.count(); }
        bool is_valid() const { return m_view.is_valid(); }
        
        // accessors
        TT const *ptr(int row) const { return m_view.ptr(row); }
        TT *ptr(int row) { return m_view.ptr(row); }
        TT const &operator()(int row) const { return m_view.at(row); }
        TT &operator()(int row) { return m_view.at(row); }
        TT const &at(int row) const { return m_view.at(row); }
        TT &at(int row) { return m_view.at(row); }
        TT const &index(int index) const { return m_view.index(index); }
        TT &index(int index) { return m_view.index(index); }
        
        // Return a deep copy of this object
        Vec copy() const
        {
            Vec<TT> v(rows());
            core_block_copy(v.view(), m_view);
            return v;
        }
        
        // Copy over data from another vector
        void copy_from(VecView<TT> const &vv)
        {
            core_block_copy(m_view, vv);
        }
        
        void copy_from(Vec const &v)
        {
            core_block_copy(m_view, v.view());
        }
        
        // Copy over the source elements indicated by index list
        void copy_from(Vec<TT> const &src, std::vector<int> const &index_list)
        {
            m_view.copy_from(src.view(),index_list);
        }
        
        // Paste over data from another view to a specific location.
        // This is done using block(...).copy_from(v.block(...)), but correctly handles mismatched
        // sizes and source/destination outside bounds, e.g. negative dstrowstart.
        
        void paste_from(Vec const &v, int dstrowstart)
        {
            m_view.paste_from(v.view(), dstrowstart);
        }
        
        // Fill with value
        void set(TT const &f) { m_view.set(f); }
        Vec &operator=(TT const &f) { m_view.set(f); return *this; }
        
        // Assign the elements indicated by index list
        void set(TT const &val, std::vector<int> const &index_list)
        {
            m_view.set(val, index_list);
        }
        
        // Create from matlab format text field, e.g. "[ 1 2 3 4 5 6 ]" or "[3+6i; 4-5i]"
        // Size is determined automatically.
        Vec &operator<<(char const *ptext)
        {
            Mtx<TT> m(ptext);
            int size = std::max(m.rows(), m.cols());
            resize(size);
            if(m.cols()>m.rows())
                m_view.copy_from(m.view().row(0));
            else
                m_view.copy_from(m.view().col(0));
            
            return *this;
        }
        
        // Returns a square matrix with the vector elements as the diagonal
        Mtx<TT> diag_matrix() const
        {
            Mtx<TT> m(rows(),rows());
            m = (TT)0;
            m.view().diag().copy_from(m_view);
            return m;
        }
        
        // Conversion of length 1 vector to TT
        operator TT() const
        {
            if(rows()!=1)
                IM_THROW_MATRIX;
            
            return at(0);
        }
        
        // Vector sub view
        Vec const block(int rowstart,
                        int rowcount,
                        int rowstep = 1) const
        {
            return Vec(m_mem, m_view.block(rowstart, rowcount, rowstep));
        }
        
        Vec block(int rowstart,
                  int rowcount,
                  int rowstep = 1)
        {
            return Vec(m_mem, m_view.block(rowstart, rowcount, rowstep));
        }
        
        // First few rows
        Vec const head(int nrows) const
        {
            return Vec(m_mem, m_view.head(nrows));
        }
        
        Vec head(int nrows)
        {
            return Vec(m_mem, m_view.head(nrows));
        }
        
        // Last few rows
        Vec const tail(int nrows) const
        {
            return Vec(m_mem, m_view.tail(nrows));
        }
        
        Vec tail(int nrows)
        {
            return Vec(m_mem, m_view.tail(nrows));
        }
        
        Mtx<TT> const matrix() const
        {
            return Mtx<TT>(m_mem, m_view.matrix_view());
        }
        
        Mtx<TT> matrix()
        {
            return Mtx<TT>(m_mem, m_view.matrix_view());
        }
        
        // Reshape the 1D matrix into a 2D matrix format
        
        // Row major
        Mtx<TT> const rm_2d_matrix(int colcount) const
        {
            return Mtx<TT>(m_mem, m_view.rm_2d_matrix_view(colcount));
        }
        
        Mtx<TT> rm_2d_matrix(int colcount)
        {
            return Mtx<TT>(m_mem, m_view.rm_2d_matrix_view(colcount));
        }
        
        // Column major
        Mtx<TT> const cm_2d_matrix(int rowcount) const
        {
            return Mtx<TT>(m_mem, m_view.cm_2d_matrix_view(rowcount));
        }
        
        Mtx<TT> cm_2d_matrix(int rowcount)
        {
            return Mtx<TT>(m_mem, m_view.cm_2d_matrix_view(rowcount));
        }

        Vec const reverse() const
        {
            return Vec(m_mem, m_view.reverse());
        }
        
        Vec reverse()
        {
            return Vec(m_mem, m_view.reverse());
        }
        
        // Copy over the source elements listed in the index list into consecutive rows
        void map(Vec const &src, std::vector<int> const &index_list)
        {
            m_view.map(src.view(), index_list);
        }
        
        // Operators
        
        Vec operator+() const
        {
            Vec v(rows());
            core_block_copy(v.view(), m_view);
            return v;
        }
        
        Vec operator-() const
        {
            Vec v(rows());
            core_block_scale(v.view(), m_view, (TT)-1);
            return v;
        }
        
        // Operations with a single value
        
        Vec operator+(TT const &val) const
        {
            Vec v(rows());
            core_block_offset(v.view(), m_view, val);
            return v;
        }
        
        Vec &operator+=(TT const &val)
        {
            core_block_offset(m_view, m_view, val);
            return *this;
        }
        
        Vec operator-(TT const &val) const
        {
            Vec v(rows());
            core_block_offset(v.view(), m_view, -val);
            return v;
        }
        
        Vec &operator-=(TT const &val)
        {
            core_block_offset(m_view, m_view, -val);
            return *this;
        }
        
        Vec operator*(TT const &val) const
        {
            Vec v(rows());
            core_block_scale(v.view(), m_view, val);
            return v;
        }
        
        Vec &operator*=(TT const &val)
        {
            core_block_scale(m_view, m_view, val);
            return *this;
        }
        
        Vec operator/(TT const &val) const
        {
            Vec v(rows());
            core_block_scale(v.view(), m_view, (TT)1/val);
            return v;
        }
        
        Vec &operator/=(TT const &val)
        {
            core_block_scale(m_view, m_view, (TT)1/val);
            return *this;
        }
        
        // Operations with a vector
        
        Vec operator+(Vec const &vin) const
        {
            Vec v(rows());
            core_block_add_elements(v.view(), m_view, vin.view());
            return v;
        }
        
        Vec &operator+=(Vec const &vin)
        {
            core_block_add_elements(m_view, m_view, vin.view());
            return *this;
        }
        
        Vec operator-(Vec const &vin) const
        {
            Vec v(rows());
            core_block_sub_elements(v.view(), m_view, vin.view());
            return v;
        }
        
        Vec &operator-=(Vec const &vin)
        {
            core_block_sub_elements(m_view, m_view, vin.view());
            return *this;
        }
        
        // Element-wise multiply
        Vec operator^(Vec const &vin) const
        {
            Vec v(rows());
            core_block_multiply_elements(v.view(), m_view, vin.view());
            return v;
        }
        
        Vec &operator^=(Vec const &vin)
        {
            core_block_multiply_elements(m_view, m_view, vin.view());
            return *this;
        }

        // Dot product (vector)
        TT dot_product(Vec const &vin) const
        {
            return core_block_reduce_multiply_add(m_view, vin.view());
        }
        
        // Outer product
        Mtx<TT> outer_product(Vec const &vin) const
        {
            Mtx<TT> m(rows(),vin.rows());
            core_block_outer_product(m.view(), m_view, vin.view());
            return m;
        }
        
        Vec cross_product(Vec const &vin) const
        {
            IM_CHECK_VALID(vin);
            IM_CHECK_ARGS(vin.rows()==3);
            IM_CHECK_ARGS(rows()==3);
            
            Vec v(3);
            v(0) = at(1) * vin(2) - at(2) * vin(1);
            v(1) = at(2) * vin(0) - at(0) * vin(2);
            v(2) = at(0) * vin(1) - at(1) * vin(0);
            return v;
        }
        
        Vec abs() const
        {
            Vec v(rows());
            core_block_abs(v.view(), m_view);
            return v;
        }
        
        Vec sqrt() const
        {
            Vec v(rows());
            core_block_sqrt(v.view(), m_view);
            return v;
        }
        
        Vec pow(TT const &exponent) const
        {
            Vec v(rows());
            core_block_pow(v.view(), m_view, exponent);
            return v;
        }
        
        Vec exp() const
        {
            Vec v(rows());
            core_block_exp(v.view(), m_view);
            return v;
        }
        
        Vec log() const
        {
            Vec v(rows());
            core_block_log(v.view(), m_view);
            return v;
        }
        
        // 1/(1+std::exp(-x))
        Vec sigm() const
        {
            Vec v(rows());
            core_block_sigm(v.view(), m_view);
            return v;
        }
        
        // 1/x
        Vec recip(TT const &epsilon = (TT)0) const
        {
            Vec v(rows());
            core_block_recip(v.view(), m_view, epsilon);
            return v;
        }
        
        Vec scale_offset(TT const &scale, TT const &offset) const
        {
            Vec v(rows());
            core_block_scale_offset(v.view(), m_view, scale, offset);
            return v;
        }
        
        Vec blend(Vec const &mat, float alpha) const
        {
            Vec v(rows());
            core_block_blend(v.view(), m_view, mat.view(), alpha);
            return v;
        }
        
        void normalize(TT const &k = (TT)0)
        {
            core_block_normalize(m_view, m_view, k);
        }
        
        TT magnitude() const
        {
            return core_block_blas_nrm2(m_view);
        }
        
        TT magnitude_squared() const
        {
            TT m = core_block_blas_nrm2(m_view);
            return m*m;
        }

        Vec conj() const
        {
            Vec v(rows());
            core_block_conj(v.view(), m_view);
            return v;
        }
        
        TT add() const
        {
            return core_block_reduce_add(m_view);
        }
        
        TT multiply() const
        {
            return core_block_reduce_multiply(m_view);
        }
        
        TT min() const
        {
            return core_block_reduce_min(m_view);
        }
        
        TT max() const
        {
            return core_block_reduce_max(m_view);
        }
        
        TT min_abs() const
        {
            return core_block_reduce_min_abs(m_view);
        }
        
        TT max_abs() const
        {
            return core_block_reduce_max_abs(m_view);
        }
        
        TT median(int k = -1) const
        {
            return core_block_reduce_median(m_view, k);
        }
        
        TT mean() const
        {
            return core_block_reduce_mean(m_view);
        }
        
        TT variance(TT *pmeanrtn = NULL) const
        {
            return core_block_reduce_variance(m_view, pmeanrtn);
        }
        
        TT distance(Vec const &vin) const
        {
            return std::sqrt(core_block_reduce_squared_distance(m_view, vin.view()));
        }
        
        TT distance_squared(Vec const &vin) const
        {
            return core_block_reduce_squared_distance(m_view, vin.view());
        }
        
        int location_of_max() const { return core_block_location_of_max(m_view); }
        int location_of_max_abs() const { return core_block_location_of_max_abs(m_view); }
        int location_of_min() const { return core_block_location_of_min(m_view); }
        int location_of_min_abs() const { return core_block_location_of_min_abs(m_view); }
        
        void sort(SortDirection direction = SortDirectionAscending)
        {
            core_block_sort(m_view, direction);
        }
        
        void sort_abs(SortDirection direction = SortDirectionAscending)
        {
            core_block_sort_abs(m_view, direction);
        }
        
        std::vector<int> sort_to_index_list(SortDirection direction = SortDirectionAscending) const
        {
            std::vector<int> index_list;
            core_block_sort(index_list, m_view, direction);
            return index_list;
        }
        
        std::vector<int> sort_to_index_list_abs(SortDirection direction = SortDirectionAscending) const
        {
            std::vector<int> index_list;
            core_block_sort_abs(index_list, m_view, direction);
            return index_list;
        }
        
        // Returns index list for values matching comparison
        std::vector<int> select_equal_to(TT const &val);
        std::vector<int> select_greater_than(TT const &val);
        std::vector<int> select_less_than(TT const &val);
        std::vector<int> select_greater_than_abs(TT const &val); // compare using std::abs
        std::vector<int> select_less_than_abs(TT const &val); // compare using std::abs
        
        // Returns index list for values matching comparison against elements of equal size matrix
        std::vector<int> select_equal_to(Vec const &src);
        std::vector<int> select_greater_than(Vec const &src);
        std::vector<int> select_less_than(Vec const &src);
        std::vector<int> select_greater_than_abs(Vec const &src); // compare using std::abs
        std::vector<int> select_less_than_abs(Vec const &src); // compare using std::abs
        
        // In the case of complex numbers, this only sets the real part
        void random_uniform(Rand &rnd, TT const &low = TT(0), TT const &high = TT(1));
        void random_gaussian(Rand &rnd, TT const &mean = TT(0), TT const &stddev = TT(1));
        
        TT sample_bilinear(float row) const;
        TT sample_bicubic(float row) const;
        
        //
        
        void print_size(FILE *fp = stdout) const
        {
            m_view.print_size(fp);
        }
        
        // Print in Matlab format if known type
        void print(FILE *fp = stdout) const
        {
            m_view.print(fp);
        }
        
        // Used to save as CSV etc.
        void write_text(FILE *fp, char delimiter = ' ') const
        {
            m_view.write_text(fp, delimiter);
        }
        
        // Save as a file using WriteText
        void save_text(char const *pfile, char delimiter = ' ') const
        {
            m_view.save_text(pfile, delimiter);
        }
        
        // These functions save and load the vector using the standard library data format
        void write_binary(FILE *fp) const
        {
            core_write_imp_file(fp, m_view);
        }
        
        void save_binary(char const *pfile) const
        {
            core_save_imp_file(pfile, m_view);
        }
        
        void read_binary(FILE *fp)
        {
            core_read_imp_file(fp, m_view);
        }
        
        void load_binary(char const *pfile)
        {
            core_load_imp_file(pfile, m_view);
        }
        
    };
 
    template <typename TT> Vec<TT> operator+(TT const &val, Vec<TT> const &vin)
    {
        Vec<TT> v(vin.rows());
        core_block_offset(v.view(), vin.view(), val);
        return v;
    }
    
    template <typename TT> Mtx<TT> operator-(TT const &val, Vec<TT> const &vin)
    {
        Vec<TT> v(vin.rows());
        core_block_scale_offset(v.view(), vin.view(), (TT)-1, val);
        return v;
    }
    
    template <typename TT> Mtx<TT> operator*(TT const &val, Vec<TT> const &vin)
    {
        Vec<TT> v(vin.rows());
        core_block_scale(v.view(), vin.view(), val);
        return v;
    }
    
    typedef Vec<float> Vf;
    typedef Vec<double> Vd;
    typedef Vec<im::Cf> Vcf;
    typedef Vec<im::Cd> Vcd;
    
}

#endif
