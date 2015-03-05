//
//  stats.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/26/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_stats_h
#define Metaphor_stats_h

namespace im
{
    // Compute a simple histogram vector from all the data in an array
    void core_histogram(std::vector<int> &hist_rtn, MVf const &mavsrc, float minval, float maxval, int bins);
    
    // Compute Principal Components of the data
    // Each data sample is expected to be a set of columns, one for each data vector and the number
    // of rows is the dimensionality
    // On output, the columns of mV are the principal vectors and vL is a vector of principal values
    // Vectors and values are sorted so that largest principal values come first. mV is orthogonal
    // This method calculates the covariance matrix and uses eigen decomposition
    // It is more accurate to use SVD, but there may be speed and memory issues for large data sets
    void core_stats_pca(MVf mV, VVf vL, const MVf &mData);
    
    // Fit a line through a 2D set of points
    // The observation matrix mcoordsv has two rows and n columns, where each column is an (x y) pair
    // The coefficients vector has the form (a b c), for the line equation ax + by + c = 0
    
    // Simple least squares line fit with error in the y direction only
    // The coordinate matrix has two or three rows. If it has three rows then the third row gives
    // the standard deviation for each point
    template <typename TT> void core_line_fit_2d(VecView<TT> vcoefs, MtxView<TT> const &mcoords);

    // Robust orthogonal least squares with outlier rejection using RANSAC.
    template <typename TT> void core_line_fit_2d_orthog_robust(VecView<TT> vcoefs, MtxView<TT> const &mcoords, Rand &rnd);

    // Orthogonal least squares
    // This finds a line in n-dimensions that passes least squares optimally throught through the given points with equal error on all axes
    // Points are expressed in mcoords as k columns of n dimensional column vectors
    // Result is expressed using the line equation x = vorigin + lambda * vunit
    template <typename TT> void core_line_fit_orthog(VecView<TT> vorigin, VecView<TT> vunit, MtxView<TT> const &mcoords);
    
    // Solves the least squares system A x = b where the matrix A has more rows than columns.
    // This is done using LDL decomposition of A^T A and back substitution with A^T b to solve the normal equations.
    // This method is quick for simple systems, but more precise results can be obtained for less well behaved systems by
    // using the SVD method because A^T A can be close to singular in some systems.
    // A is the design matrix where each row is one observation [X0(xi)/si .. Xn(xi)/si], where si are the
    // standard deviations (if used), and b is the vector [ yi/si ].
    // The model is given by y(x) = sum_k { a_k X_k(x) }, and a_k are the coefficients to be determined (vx).
    template <typename TT> void core_least_squares_normal(VecView<TT> vx, const MtxView<TT> &mA, const VecView<TT> &vb);
    
    // Solves the least squares system using SVD, slower but more robust when data are close to being degenerate.
    template <typename TT> void core_least_squares_svd(VecView<TT> vx, const MtxView<TT> &mA, const VecView<TT> &vb);
    
    
    // This class estimates the mean vector and full covariance matrix of multidimensional data
    // Vectors can be added sequantially and the covariance calculated at any time
    class CovarianceEstimator
    {
    public:
        CovarianceEstimator() : m_count(0) {}
        CovarianceEstimator(MVf const &mavdata) { start(mavdata.cols()); add(mavdata); }
        
        // Resets the statistics and sets up for the designated dimensionality
        void start(int data_dimension)
        {
            IM_CHECK_LOWER_BOUNDS(data_dimension, 1);
            
            m_matsumcov.resize(data_dimension, data_dimension);
            m_matsummean.resize(data_dimension);
            
            m_matsumcov = 0.0f;
            m_matsummean = 0.0f;
            m_count = 0;
        }
        
        // Add data vectors in column format
        // Each data vector is a column of mtx, and the number of rows must equal the dimensionality
        // Use matrix transpose to add data in row format
        void add(MVf const &m);
        void add(VVf const &v);
        
        // Returns the current covariance matrix
        Mf covariance() const;
        
        // Returns the current correlation matrix
        Mf correlation() const;
        
        // Returns the current mean vector
        Vf mean() const;
        
        // Returns the count of data vectors that were added
        int count() const { return m_count; }
        int dimension() const { return m_matsummean.rows(); }
        
    private:
        Md m_matsumcov;
        Vd m_matsummean;
        int m_count;
    };
    
    // This class is used to create a set of test vectors that have Gaussian statistics and are distributed
    // according to the given mean and covariance parameters.
    class GaussianNoiseVector
    {
    public:
        GaussianNoiseVector() {}
        GaussianNoiseVector(VVf const &vmean, MVf const &mcov) { init(vmean, mcov); }
        
        void init(VVf const &vmean, MVf const &mcov);
        
        // Gets a set of rows, where each column is a data vector sample
        // Pass in a random number generator
        void create_data(Rand &rnd, MVf mdata);
        
        int dimension() { return (int)m_vmean.rows(); }
        
    private:
        Mf m_matA;
        Vf m_vmean;
    };

    // Gaussian in N dimensions
    class GaussianSpace
    {
    public:
        GaussianSpace(VVf const &vmean, MVf const &mcov)
        {
            init(vmean, mcov);
        }
        
        // init with mean and covariance of Gaussian
        void init(VVf const &vmean, MVf const &mcov);
        
        // init with mean and variance of Gaussian
        void init(VVf const &vmean, float variance);
        
        // sample the Gaussian function
        double get(const VVf &vx);
        
        // compute the Mahalanobis distance between two points
        double distance(const VVf &v1,  const VVf &v2);
        
    private:
        int m_size;
        Vf m_vmean;
        MatrixDecompLDLT<float> m_ldlt;
        double m_const;
        float m_variance;
        bool m_use_variance;
    };
    
    
}

#endif
