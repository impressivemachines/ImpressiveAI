//
//  matrix.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/29/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_matrix_h
#define Metaphor_matrix_h

namespace im
{
    template <typename TT> class Vec;
    
    // While MtxView is just a view on memory, Mtx also manages memory.
    // Mtx also defines a number of additional operations such as matrix operators +* etc.
    // Mtx is only implemented for float, double, complex float and complex double.
    // Mtx can wrap external memory by initializing it from a MtxView.
    
    template <typename TT>
    class Mtx
    {
    private:
        MtxView<TT> m_view;
        std::shared_ptr<std::vector<TT>> m_mem;
      
    protected:
        friend class Vec<TT>;
        
        Mtx(std::shared_ptr<std::vector<TT>> const &mem, MtxView<TT> const &mav) : m_mem(mem), m_view(mav) {}
        
    public:
        typedef TT ValueType;
        
        Mtx() { resize(0, 0); }
        
        // Construct a matrix
        Mtx(int nrows, int ncols) { resize(nrows, ncols); }
        
        // Construct using << operator from matlab formatted text string representation, eg. "[ 1 2; 3 4 ]"
        Mtx(char const *pmatlabtext) { *this << pmatlabtext; }
        
        // Wrap the given view - note that this object then manages no memory and that there is no data copy
        Mtx(MtxView<TT> const &mav) { m_mem.reset(); m_view = mav; }
        
        Mtx(int rows, int cols, int rowstride, int colstride, TT const *pdata)
        {
            m_mem.reset();
            m_view.wrap(rows, cols, rowstride, colstride, pdata);
        }
        
        TT const *wrap(int nrows, int ncols, int rowstride, int colstride, TT const *pdata)
        {
            m_mem.reset();
            return m_view.wrap(nrows, ncols, rowstride, colstride, pdata);
        }
        
        void wrap(const std::vector<TT> &v)
        {
            m_mem.reset();
            m_view.wrap(v);
        }
        
        // Resizes the matrix and loses any data
        void resize(int nrows, int ncols);
        
        // Wrap the given view - note that this object then manages no memory and that there is no data copy
        Mtx &operator=(MtxView<TT> const &mav) { m_mem.reset(); m_view = mav; return *this; }
        
        // Return the view
        MtxView<TT> const &view() const { return m_view; }
        MtxView<TT> view() { return m_view; }
        
        void deallocate() { m_view.reset(); m_mem.reset(); }
        
        // Stops sharing the matrix by creating own private copy and guarentees packed row major layout
        void stop_sharing();
        
        // Ensures matrix has the requested capacity non-destructively
        // This allows future resize operations to proceed without allocation overhead
        void reserve(int nrows, int ncols)
        {
            IM_CHECK_LOWER_BOUNDS(nrows, 0);
            IM_CHECK_LOWER_BOUNDS(ncols, 0);
            
            if(m_mem)
                m_mem->reserve(nrows * ncols);
        }
        
        // Add one or more rows at the bottom of the matrix
        void push_back(MtxView<TT> const &mav);
        
        // Add an new row at the bottom of the matrix and sets its element(s) to the given value
        void push_back(TT const &f);
        
        // Remove a row from the bottom of the matrix
        void pop_back();
        
        // Access information
        int rows() const { return m_view.rows(); }
        int cols() const { return m_view.cols(); }
        int row_stride() const { return m_view.row_stride(); }
        int col_stride() const { return m_view.col_stride(); }
        int count() const { return m_view.count(); }
        bool is_valid() const { return m_view.is_valid(); }
        
        // 2D accessors
        TT const *ptr(int row, int col) const { return m_view.ptr(row, col); }
        TT *ptr(int row, int col) { return m_view.ptr(row, col); }
        TT const &operator()(int row, int col) const { return m_view.at(row, col); }
        TT &operator()(int row, int col) { return m_view.at(row, col); }
        TT const &at(int row, int col) const { return m_view.at(row, col); }
        TT &at(int row, int col) { return m_view.at(row, col); }
        
        // Index as linear array regardless of shape
        TT const &index(int index) const { return m_view.index(index); }
        TT &index(int index) { return m_view.index(index); }
        
        // Return a deep copy of this object
        Mtx copy() const
        {
            Mtx m(rows(),cols());
            core_block_copy(m.view(), m_view);
            return m;
        }
        
        // Copy over data from another matrix
        void copy_from(MtxView<TT> const &mv)
        {
            core_block_copy(m_view, mv);
        }
        
        void copy_from(Mtx const &mat)
        {
            core_block_copy(m_view, mat.view());
        }
        
        // Copy only the elements from the source indicated by index list
        void copy_from(Mtx const &src, std::vector<MtxLoc> const &index_list)
        {
            copy_from(src.view(), index_list);
        }
        
        // Copy data in raster order without regard to the relative shapes.
        // Any destination elements which have no source data are set to zero.
        void reshape_from(Mtx const &mat)
        {
            if(rows()>0 && cols()>0)
                core_block_reshape(m_view, mat.view());
        }
        
        // Paste over data from another view to a specific location.
        // This is done using block(...).copy_from(mav.block(...)), but correctly handles mismatched
        // sizes and source/destination outside bounds, e.g. negative dstrowstart.
        void paste_from(Mtx const &mat, int dstrowstart, int dstcolstart)
        {
            m_view.paste_from(mat.view(), dstrowstart, dstcolstart);
        }
        
        // Fill with value
        void set(TT const &f) { m_view.set(f); }
        Mtx &operator=(TT const &f) { m_view.set(f); return *this; }
        
        // Assign the value to the elements indicated by index list
        void set(TT const &val, std::vector<MtxLoc> const &index_list)
        {
            m_view.set(val, index_list);
        }
        
        // Create from matlab format text field, e.g. "[ 1 2 3; 4 5 6 ]" or "[3+6i; 4-5i]"
        // Size is determined automatically.
        Mtx &operator<<(char const *ptext);
        
        // Identity
        void set_identity() { m_view = (TT)0; m_view.diag() = (TT)1; }
        
        // Conversion of 1x1 matrix to TT
        operator TT() const
        {
            if(rows()!=1 || cols()!=1)
                IM_THROW_MATRIX;
            
            return at(0,0);
        }
        
        // Matrix sub view
        Mtx const block(int rowstart, int colstart,
                        int rowcount, int colcount,
                        int rowstep = 1, int colstep = 1) const
        {
            return Mtx(m_mem, m_view.block(rowstart, colstart, rowcount, colcount, rowstep, colstep));
        }
        
        Mtx block(int rowstart, int colstart,
                  int rowcount, int colcount,
                  int rowstep = 1, int colstep = 1)
        {
            return Mtx(m_mem, m_view.block(rowstart, colstart, rowcount, colcount, rowstep, colstep));
        }
        
        Mtx const block(MtxRect const &rct) const { return block(rct.origin.row, rct.origin.col, rct.size.rows, rct.size.cols); }
        Mtx block(MtxRect const &rct) { return block(rct.origin.row, rct.origin.col, rct.size.rows, rct.size.cols); }
        
        // Return a row as a vector
        Vec<TT> const row(int r) const
        {
            return Vec<TT>(m_mem, m_view.row(r));
        }
        
        Vec<TT> row(int r)
        {
            return Vec<TT>(m_mem, m_view.row(r));
        }
        
        // Return a column as a vector
        Vec<TT> const col(int c) const
        {
            return Vec<TT>(m_mem, m_view.col(c));
        }
        
        Vec<TT> col(int c)
        {
            return Vec<TT>(m_mem, m_view.col(c));
        }
        
        // Return the diagonal as a column vector
        // Positive offset gets the diagonal starting from (0,offset)
        // Negative offset gets the diagonal starting from (-offset,0)
        Vec<TT> const diag(int offset = 0) const
        {
            return Vec<TT>(m_mem, m_view.diag(offset));
        }
        
        Vec<TT> diag(int offset = 0)
        {
            return Vec<TT>(m_mem, m_view.diag(offset));
        }
        
        // Return the transpose matrix or vector
        Mtx const t() const
        {
            return Mtx(m_mem, m_view.t());
        }
        
        Mtx t()
        {
            return Mtx(m_mem, m_view.t());
        }
        
        // Return this matrix flipped horizontally
        Mtx const reverse_cols() const
        {
            return Mtx(m_mem, m_view.reverse_cols());
        }
        
        Mtx reverse_cols()
        {
            return Mtx(m_mem, m_view.reverse_cols());
        }
        
        // Return this matrix flipped vertically
        Mtx const reverse_rows() const
        {
            return Mtx(m_mem, m_view.reverse_rows());
        }
        
        Mtx reverse_rows()
        {
            return Mtx(m_mem, m_view.reverse_rows());
        }
        
        Mtx const rotate_90cw() const
        {
            return Mtx(m_mem, m_view.rotate_90cw());
        }
        
        Mtx rotate_90cw()
        {
            return Mtx(m_mem, m_view.rotate_90cw());
        }
        
        Mtx const rotate_90ccw() const
        {
            return Mtx(m_mem, m_view.rotate_90ccw());
        }
        
        Mtx rotate_90ccw()
        {
            return Mtx(m_mem, m_view.rotate_90ccw());
        }
        
        Mtx const rotate_180() const
        {
            return Mtx(m_mem, m_view.rotate_180());
        }
        
        Mtx rotate_180()
        {
            return Mtx(m_mem, m_view.rotate_180());
        }
        
        // Copy over the source elements listed in the index list into consecutive cols/rows of *this.
        void map(Mtx const &src, std::vector<MtxLoc> const &index_list)
        {
            m_view.map(src.view(), index_list);
        }
        
        // Copy over the source rows listed in the index list into consecutive rows of *this.
        void map_rows(Mtx const &src, std::vector<int> const &index_list)
        {
            m_view.map_rows(src.view(), index_list);
        }
        
        // Copy over the source cols listed in the index list into consecutive cols of *this.
        void map_cols(Mtx const &src, std::vector<int> const &index_list)
        {
            m_view.map_cols(src.view(), index_list);
        }
        
        // Unary ops
        
        Mtx operator+() const
        {
            Mtx matrtn(rows(), cols());
            core_block_copy(matrtn.view(), m_view);
            return matrtn;
        }
        
        Mtx operator-() const
        {
            Mtx matrtn(rows(), cols());
            core_block_scale(matrtn.view(), m_view, (TT)-1);
            return matrtn;
        }
        
        // Operations with a single value
        
        Mtx operator+(TT const &val) const
        {
            Mtx matrtn(rows(), cols());
            core_block_offset(matrtn.view(), m_view, val);
            return matrtn;
        }
        
        Mtx &operator+=(TT const &val)
        {
            core_block_offset(m_view, m_view, val);
            return *this;
        }
        
        Mtx operator-(TT const &val) const
        {
            Mtx matrtn(rows(), cols());
            core_block_offset(matrtn.view(), m_view, -val);
            return matrtn;
        }
        
        Mtx &operator-=(TT const &val)
        {
            core_block_offset(m_view, m_view, -val);
            return *this;
        }
        
        Mtx operator*(TT const &val) const
        {
            Mtx matrtn(rows(), cols());
            core_block_scale(matrtn.view(), m_view, val);
            return matrtn;
        }
        
        Mtx &operator*=(TT const &val)
        {
            core_block_scale(m_view, m_view, val);
            return *this;
        }
        
        Mtx operator/(TT const &val) const
        {
            Mtx matrtn(rows(), cols());
            core_block_scale(matrtn.view(), m_view, (TT)1/val);
            return matrtn;
        }
        
        Mtx &operator/=(TT const &val)
        {
            core_block_scale(m_view, m_view, (TT)1/val);
            return *this;
        }
        
        // Operations with a matrix
        
        Mtx operator+(Mtx const &mat) const
        {
            Mtx matrtn(rows(), cols());
            core_block_add_elements(matrtn.view(), m_view, mat.view());
            return matrtn;
        }
        
        Mtx &operator+=(Mtx const &mat)
        {
            core_block_add_elements(m_view, m_view, mat.view());
            return *this;
        }
        
        Mtx operator-(Mtx const &mat) const
        {
            Mtx matrtn(rows(), cols());
            core_block_sub_elements(matrtn.view(), m_view, mat.view());
            return matrtn;
        }
        
        Mtx &operator-=(Mtx const &mat)
        {
            core_block_sub_elements(m_view, m_view, mat.view());
            return *this;
        }
        
        Mtx operator*(Mtx const &mat) const
        {
            Mtx matrtn(rows(),mat.cols());
            core_block_blas_gemm(matrtn.view(), m_view, mat.view(), (TT)1, (TT)0, TransMode_N, TransMode_N);
            return matrtn;
        }
        
        Mtx &operator*=(Mtx const &mat)
        {
            Mtx matrtn(rows(),mat.cols());
            core_block_blas_gemm(matrtn.view(), m_view, mat.view(), (TT)1, (TT)0, TransMode_N, TransMode_N);
            core_block_copy(m_view, matrtn.view());
            return *this;
        }
        
        Vec<TT> operator*(Vec<TT> const &v)
        {
            IM_CHECK_ARGS(v.rows()==cols());
            Vec<TT> vrtn(rows());
            core_block_blas_gemv(vrtn.view(), m_view, v.view(), (TT)1, (TT)0, TransMode_N);
            return vrtn;
        }
        
        // Element-wise multiply
        Mtx operator^(Mtx const &mat) const
        {
            Mtx matrtn(rows(),cols());
            core_block_multiply_elements(matrtn.view(), m_view, mat.view());
            return matrtn;
        }
        
        Mtx &operator^=(Mtx const &mat)
        {
            core_block_multiply_elements(m_view, m_view, mat.view());
            return *this;
        }
        
        // Dot product (matrix)
        TT dot_product(Mtx const &mat) const
        {
            return core_block_reduce_multiply_add(m_view, mat.view());
        }
        
        // Compute xT A x, where x is a vector and A is square
        TT xTAx(Vec<TT> const &x) const
        {
            return core_block_compute_xTAx(m_view, x.view());
        }
        
        Mtx abs() const
        {
            Mtx matrtn(rows(),cols());
            core_block_abs(matrtn.view(), m_view);
            return matrtn;
        }
        
        Mtx sqrt() const
        {
            Mtx matrtn(rows(),cols());
            core_block_sqrt(matrtn.view(), m_view);
            return matrtn;
        }
        
        Mtx pow(TT const &exponent) const
        {
            Mtx matrtn(rows(),cols());
            core_block_pow(matrtn.view(), m_view, exponent);
            return matrtn;
        }
        
        Mtx exp() const
        {
            Mtx matrtn(rows(),cols());
            core_block_exp(matrtn.view(), m_view);
            return matrtn;
        }
        
        Mtx log() const
        {
            Mtx matrtn(rows(),cols());
            core_block_log(matrtn.view(), m_view);
            return matrtn;
        }
        
        // 1/(1+std::exp(-x))
        Mtx sigm() const
        {
            Mtx matrtn(rows(),cols());
            core_block_sigm(matrtn.view(), m_view);
            return matrtn;
        }
        
        // 1/x
        Mtx recip(TT const &epsilon = (TT)0) const
        {
            Mtx matrtn(rows(),cols());
            core_block_recip(matrtn.view(), m_view, epsilon);
            return matrtn;
        }
        
        Mtx scale_offset(TT const &scale, TT const &offset) const
        {
            Mtx matrtn(rows(),cols());
            core_block_scale_offset(matrtn.view(), m_view, scale, offset);
            return matrtn;
        }
        
        Mtx blend(Mtx const &mat, float alpha) const
        {
            Mtx matrtn(rows(),cols());
            core_block_blend(matrtn.view(), m_view, mat.view(), alpha);
            return matrtn;
        }
        
        void normalize(TT const &k = (TT)0)
        {
            core_block_normalize(m_view, m_view, k);
        }
        
        void normalize_rows(TT const &k = (TT)0)
        {
            MtxView<TT> v = m_view.t(); // transpose
            core_block_normalize_columns(v, v, k);
        }
        
        void normalize_columns(TT const &k = (TT)0)
        {
            core_block_normalize_columns(m_view, m_view, k);
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
        
        TT trace() const
        {
            return core_block_reduce_add(m_view.diag());
        }
        
        Mtx conj() const
        {
            Mtx matrtn(rows(),cols());
            core_block_conj(matrtn.view(), m_view);
            return matrtn;
        }
        
        Mtx conj_t() const
        {
            Mtx matrtn(cols(),rows());
            core_block_conj(matrtn.view(), m_view.t());
            return matrtn;
        }
        
        Mtx inverse() const;
        
        TT det() const;
        
        TT log_det() const;
        
        TT sign_det() const;
        
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
        
        TT distance(Mtx const &mat) const
        {
            return std::sqrt(core_block_reduce_squared_distance(m_view, mat.view()));
        }
        
        TT distance_squared(Mtx const &mat) const
        {
            return core_block_reduce_squared_distance(m_view, mat.view());
        }
        
        // Reductions along each row - use t() to reduce columns
        
        Vec<TT> reduce_rows_add() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_add(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_multiply() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_multiply(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_min() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_min(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_max() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_max(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_min_abs() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_min_abs(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_max_abs() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_max_abs(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_median(int k = -1) const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_median(v.view(), m_view, k);
            return v;
        }
        
        Vec<TT> reduce_rows_mean() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_mean(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_variance() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_variance(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_sum_squares() const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_sum_squares(v.view(), m_view);
            return v;
        }
        
        Vec<TT> reduce_rows_distance_squared(Mtx const &mat) const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_squared_distance(v.view(), m_view, mat.view());
            return v;
        }
        
        Vec<TT> reduce_rows_dot_product(Mtx const &mat) const
        {
            Vec<TT> v(rows());
            core_block_reduce_rows_multiply_add(v.view(), m_view, mat.view());
            return v;
        }
        
        MtxLoc location_of_max() const { return core_block_location_of_max(m_view); }
        MtxLoc location_of_max_abs() const { return core_block_location_of_max_abs(m_view); }
        MtxLoc location_of_min() const { return core_block_location_of_min(m_view); }
        MtxLoc location_of_min_abs() const { return core_block_location_of_min_abs(m_view); }
        
        void sort(SortDirection direction = SortDirectionAscending)
        {
            core_block_sort(m_view, direction);
        }
        
        void sort_abs(SortDirection direction = SortDirectionAscending)
        {
            core_block_sort_abs(m_view, direction);
        }
        
        std::vector<MtxLoc> sort_to_index_list(SortDirection direction = SortDirectionAscending) const
        {
            std::vector<MtxLoc> index_list;
            core_block_sort(index_list, m_view, direction);
            return index_list;
        }
        
        std::vector<MtxLoc> sort_to_index_list_abs(SortDirection direction = SortDirectionAscending) const
        {
            std::vector<MtxLoc> index_list;
            core_block_sort_abs(index_list, m_view, direction);
            return index_list;
        }

        // Returns index list for values matching comparison
        std::vector<MtxLoc> select_equal_to(TT const &val);
        std::vector<MtxLoc> select_greater_than(TT const &val);
        std::vector<MtxLoc> select_less_than(TT const &val);
        std::vector<MtxLoc> select_greater_than_abs(TT const &val); // compare using std::abs
        std::vector<MtxLoc> select_less_than_abs(TT const &val); // compare using std::abs
        
        // Returns index list for values matching comparison against elements of equal size matrix
        std::vector<MtxLoc> select_equal_to(Mtx const &src);
        std::vector<MtxLoc> select_greater_than(Mtx const &src);
        std::vector<MtxLoc> select_less_than(Mtx const &src);
        std::vector<MtxLoc> select_greater_than_abs(Mtx const &src); // compare using std::abs
        std::vector<MtxLoc> select_less_than_abs(Mtx const &src); // compare using std::abs
        
        // In the case of complex numbers, this only sets the real part
        void random_uniform(Rand &rnd, TT const &low = TT(0), TT const &high = TT(1));
        void random_gaussian(Rand &rnd, TT const &mean = TT(0), TT const &stddev = TT(1));
        
        TT sample_bilinear(float row, float col) const;
        TT sample_bicubic(float row, float col) const;

        //
        
        void print_size(FILE *fp = stdout) const
        {
            m_view.print_size(fp);
        }
        
        // Print in Matlab format
        void print(FILE *fp = stdout) const
        {
            m_view.print(fp);
        }
        
        // Writes just the data with no size information as a text file, one row per line
        // In this format, complex numbers are split into separate real and imaginary fields
        // Use delimiter = ',' for CSV format
        void write_text(FILE *fp, char delimiter = ' ') const
        {
            m_view.write_text(fp, delimiter);
        }
        
        void save_text(char const *pfile, char delimiter = ' ') const
        {
            m_view.save_text(pfile, delimiter);
        }
        
        // These functions save and load the matrix using the standard library data format
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
    
    template <typename TT> Mtx<TT> operator+(TT const &val, Mtx<TT> const &mat)
    {
        Mtx<TT> matrtn(mat.rows(), mat.cols());
        core_block_offset(matrtn.view(), mat.view(), val);
        return matrtn;
    }
    
    template <typename TT> Mtx<TT> operator-(TT const &val, Mtx<TT> const &mat)
    {
        Mtx<TT> matrtn(mat.rows(), mat.cols());
        core_block_scale_offset(matrtn.view(), mat.view(), (TT)-1, val);
        return matrtn;
    }
    
    template <typename TT> Mtx<TT> operator*(TT const &val, Mtx<TT> const &mat)
    {
        Mtx<TT> matrtn(mat.rows(), mat.cols());
        core_block_scale(matrtn.view(), mat.view(), val);
        return matrtn;
    }
    
    
    typedef Mtx<float> Mf;
    typedef Mtx<double> Md;
    typedef Mtx<im::Cf> Mcf;
    typedef Mtx<im::Cd> Mcd;
}

#endif
