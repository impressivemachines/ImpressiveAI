//
//  stats.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/26/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

void im::core_histogram(std::vector<int> &hist_rtn, MVf const &mavsrc, float minval, float maxval, int bins)
{
    IM_CHECK_LOWER_BOUNDS(bins, 1);
    IM_CHECK_ARGS(maxval > minval);

    hist_rtn.resize(bins);
    
    memset(hist_rtn.data(), 0, bins * sizeof(int));
    
    float scale = (bins - 1)/(maxval - minval);
    
    int r,c;
    for(r=0; r<mavsrc.rows(); r++)
        for(c=0; c<mavsrc.cols(); c++)
        {
            int ind = core_round(scale * (mavsrc(r,c) - minval));
            if(ind<0)
                ind = 0;
            else if(ind>=bins)
                ind = bins-1;
            hist_rtn[ind]++;
        }
}

/*template <typename TT>
bool im::core_least_squares_2d(VecView<TT> vcoefs, MtxView<TT> const &mcoords)
{
    return true;
}

template <typename TT>
bool im::core_least_squares_normal(VecView<TT> vx, const MtxView<TT> &mA, const VecView<TT> &vb)
{
    return true;
}

template <typename TT>
bool im::core_least_squares_svd(VecView<TT> vx, const MtxView<TT> &mA, const VecView<TT> &vb)
{
    return true;
}
*/

void im::CovarianceEstimator::add(MVf const &m)
{
    int d = dimension();
    IM_CHECK_VALID(m);
    IM_CHECK_ARGS(m.rows()==d);
    
    for(int col=0; col<m.cols(); col++)
        for(int i=0; i<d; i++)
        {
            for(int j=0; j<=i; j++)
                m_matsumcov(i,j) += (double)m(i,col) * (double)m(j,col);
            m_matsummean(i) += m(i,col);
        }
    
    m_count+=m.cols();
}

void im::CovarianceEstimator::add(VVf const &v)
{
    int d = dimension();
    IM_CHECK_VALID(v);
    IM_CHECK_ARGS(v.rows()==d);
    
    for(int i=0; i<d; i++)
    {
        for(int j=0; j<=i; j++)
            m_matsumcov(i,j) += (double)v(i) * (double)v(j);
        m_matsummean(i) += v(i);
    }
    
    m_count++;
}

// Returns the current covariance matrix
im::Mf im::CovarianceEstimator::covariance() const
{
    IM_CHECK_LOWER_BOUNDS(m_count, 2);
    int d = dimension();
    Mf mcov(d,d);
    
    for(int i=0; i<d; i++)
        for(int j=0; j<=i; j++)
            mcov(i,j) = (float)(m_matsumcov(i,j)/m_count - (m_matsummean(i)/m_count) * (m_matsummean(j)/m_count));
    
    core_block_copy_upper_tri(mcov.view(), mcov.view().t(), false);
    
    return mcov;
}

// Returns the current correlation matrix
im::Mf im::CovarianceEstimator::correlation() const
{
    IM_CHECK_LOWER_BOUNDS(m_count, 2);
    int d = dimension();
    Mf mcorr(d,d);
    
    for(int i=0; i<d; i++)
        for(int j=0; j<=i; j++)
        {
            mcorr(i,j) = (float)(m_matsumcov(i,j)/m_count - (m_matsummean(i)/m_count) * (m_matsummean(j)/m_count));
            if(i==j)
                mcorr(i,i) = sqrtf(mcorr(i,i));
        }
    
    for(int i=1; i<d; i++)
        for(int j=0; j<i; j++)
            mcorr(i,j) /= mcorr(i,i)*mcorr(j,j);
    
    for(int i=0; i<d; i++)
        mcorr(i,i) = 1;
    
    core_block_copy_upper_tri(mcorr.view(), mcorr.view().t(), false);
    
    return mcorr;
}

// Returns the current mean vector
im::Vf im::CovarianceEstimator::mean() const
{
    IM_CHECK_LOWER_BOUNDS(m_count, 1);
    int d = dimension();
    Vf vmean(d);
    
    for(int i=0; i<d; i++)
        vmean(i) = (float)(m_matsummean(i)/m_count);
    
    return vmean;
}

void im::core_stats_pca(MVf mV, VVf vL, const MVf &mData)
{
    IM_CHECK_VALID(mV);
    IM_CHECK_VALID(vL);
    IM_CHECK_VALID(mData);
    int d = mData.rows();
    IM_CHECK_MATRIX_SQUARE(mV);
    IM_CHECK_ARGS(mV.rows()==d);
    IM_CHECK_ARGS(vL.rows()==d);
    
    CovarianceEstimator cov;
    cov.start(d);
    cov.add(mData);
    
    Mf mcov = cov.covariance();
    
    MatrixDecompEigenSymmetric<float> eig(mcov.view());
    
    mV.copy_from(eig.eigenvectors().view());
    vL.copy_from(eig.eigenvalues().view());
}

void im::GaussianNoiseVector::init(VVf const &vmean, MVf const &mcov)
{
    IM_CHECK_VALID(vmean);
    IM_CHECK_VALID(mcov);
    IM_CHECK_LOWER_BOUNDS(vmean.rows(), 1);
    IM_CHECK_MATRIX_SQUARE(mcov);
    IM_CHECK_ARGS(vmean.rows()==mcov.rows());
    
    m_vmean.copy_from(vmean);
    MatrixDecompLLT<float> lltCov(mcov); // compute the Cholesky decomposition of covariance
    m_matA = lltCov.matrix_L();
}

// Gets a set of rows, where each column is a data vector sample
// Pass in a random number generator
void im::GaussianNoiseVector::create_data(Rand &rnd, MVf mdata)
{
    int d = dimension();
    IM_CHECK_VALID(mdata);
    IM_CHECK_ARGS(mdata.rows()==d);
    
    for(int i=0; i<mdata.cols(); i++)
    {
        Vf v(d);
        for(int j=0; j<d; j++)
            v(j) = rnd.gauss();
        
        v = m_matA * v + m_vmean;
        
        for(int j=0; j<d; j++)
            mdata(j,i) = v(j);
    }
}

void im::GaussianSpace::init(VVf const &vmean, MVf const &mcov)
{
    IM_CHECK_VALID(vmean);
    IM_CHECK_VALID(mcov);
    IM_CHECK_ARGS(vmean.rows()==mcov.rows());
    IM_CHECK_MATRIX_SQUARE(mcov);
    
    m_size = vmean.rows();
    m_vmean.resize(m_size);
    m_vmean.copy_from(vmean);
    m_ldlt.compute(mcov);
    m_const = 1.0/std::sqrt(std::pow(2*CONST_PI, (double)m_size) * m_ldlt.vector_D().multiply());
    m_use_variance = false;
}

void im::GaussianSpace::init(VVf const &vmean, float variance)
{
    IM_CHECK_VALID(vmean);
    IM_CHECK_ARGS(variance>0);
    
    m_size = vmean.rows();
    m_vmean.resize(m_size);
    m_vmean.copy_from(vmean);
    m_variance = variance;
    m_const = 1.0/std::sqrt(std::pow(2*CONST_PI, (double)m_size) * pow((double)variance, (double)m_size));
    m_use_variance = true;
}

double im::GaussianSpace::get(const VVf &vx)
{
    IM_CHECK_VALID(vx);
    IM_CHECK_VECTOR_SIZE(vx, m_size);
    
    Vf vdx = Vf(vx) - m_vmean;
    if(m_use_variance)
        return m_const * std::exp(-0.5 * vdx.magnitude_squared() / m_variance);
    else
        return m_const * std::exp(-0.5 * vdx.dot_product(m_ldlt.solve(vdx.view())));
}

double im::GaussianSpace::distance(const VVf &v1,  const VVf &v2)
{
    IM_CHECK_VALID(v1);
    IM_CHECK_VECTOR_SIZE(v1, m_size);
    IM_CHECK_VALID(v2);
    IM_CHECK_VECTOR_SIZE(v2, m_size);
    
    Vf vdx = Vf(v1) - Vf(v2);
    if(m_use_variance)
        return vdx.magnitude();
    else
        return std::sqrt(vdx.dot_product(m_ldlt.solve(vdx.view())));
}
