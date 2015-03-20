//
//  decomp_eigen.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT>
static void core_priv_make_givens(TT const &a, TT const &b, TT &c, TT &s)
{
    if(b==0)
    {
        c = (TT)1;
        s = (TT)0;
    }
    else if(std::abs(b) > std::abs(a))
    {
        TT t = -a / b;
        TT s1 = (TT)1.0 / std::sqrt(1 + t * t);
        s = s1;
        c = s1 * t;
    }
    else
    {
        TT t = -b / a;
        TT c1 = (TT)1.0 / std::sqrt(1 + t * t);
        c = c1;
        s = c1 * t;
    }
}

template <typename TT>
static void core_priv_chop_small(im::VecView<TT> subdiag, im::VecView<TT> const &diag)
{
    int N = diag.rows();
    for(int i=0; i<N-1; i++)
        if(std::abs(subdiag(i)) < im::TypeProperties<TT>::epsilon() * (std::abs(diag(i)) + std::abs(diag(i+1))))
            subdiag(i) = (TT)0;
}

template <typename TT>
static void core_priv_qrstep(im::VecView<TT> vdiag, im::VecView<TT> vsubdiag, im::VecView<TT> vgc, im::VecView<TT> vgs)
{
    int n = vdiag.rows();
    
    IM_CHECK_LOWER_BOUNDS(n,2);
    
    // trailing eignevalue
    TT ta = vdiag(n-2);
    TT tb = vdiag(n-1);
    TT tab = vsubdiag(n-2);
    
    TT dt = (ta - tb)/(TT)2;
    TT mu;
    
    if(dt>(TT)0)
        mu = tb - tab * (tab / (dt + im::core_hypot(dt,tab)));
    else if(dt==(TT)0)
        mu = tb - std::abs(tab);
    else
        mu = tb + tab * (tab / ((-dt) + im::core_hypot(dt,tab)));
    
    if(im::TypeProperties<TT>::epsilon() * std::abs(mu) > (std::abs(vdiag(0)) + std::abs(vsubdiag(0))))
        mu = (TT)0;
    
    TT x = vdiag(0) - mu;
    TT z = vsubdiag(0);
    TT ak = 0;
    TT bk = 0;
    TT zk = 0;
    TT ap = vdiag(0);
    TT bp = vsubdiag(0);
    TT aq = vdiag(1);
    
    if(n==2)
    {
        TT c, s;
        core_priv_make_givens(x, z, c, s);
        
        vgc(0) = c;
        vgs(0) = s;
        TT ap1 = c * (c * ap - s * bp) + s * (s * aq - c * bp);
        TT bp1 = c * (s * ap + c * bp) - s * (s * bp + c * aq);
        TT aq1 = s * (s * ap + c * bp) + c * (s * bp + c * aq);
        ak = ap1;
        bk = bp1;
        ap = aq1;
        vdiag(0) = ak;
        vsubdiag(0) = bk;
        vdiag(1) = ap;
        return;
    }
    
    TT bq = vsubdiag(1);
    
    for(int k = 0; k < n-1; k++)
    {
        TT c, s;
        core_priv_make_givens(x, z, c, s);
        
        // store Givens rotation
        vgc(k) = c;
        vgs(k) = s;
        
        // compute G' T G
        TT bk1 = c * bk - s * zk;
        TT ap1 = c * (c * ap - s * bp) + s * (s * aq - c * bp);
        TT bp1 = c * (s * ap + c * bp) - s * (s * bp + c * aq);
        TT zp1 = -s * bq;
        TT aq1 = s * (s * ap + c * bp) + c * (s * bp + c * aq);
        TT bq1 = c * bq;
        
        ak = ap1;
        bk = bp1;
        zk = zp1;
        ap = aq1;
        bp = bq1;
        
        if(k < n - 2)
            aq = vdiag(k + 2);
        if(k < n - 3)
            bq = vsubdiag(k + 2);
        
        vdiag(k) = ak;
        
        if(k > 0)
            vsubdiag(k - 1) = bk1;
        if(k < n - 2)
            vsubdiag(k + 1) = bp;
        
        x = bk;
        z = zk;
    }
    
    vdiag(n-1) = ap;
    vsubdiag(n-2) = bk;
}

template <typename TT>
void im::MatrixDecompEigenSymmetric<TT>::compute(MtxView<TT> const &mavA, bool compute_vectors)
{
    IM_CHECK_VALID(mavA);
    IM_CHECK_MATRIX_SQUARE(mavA);
    
    int size = mavA.rows();
    m_vdiag.resize(size);
    m_vsubdiag.resize(size-1);
    
    if(size==1)
    {
        m_vdiag(0) = mavA(0,0);
        m_mvectors.resize(1,1);
        m_mvectors(0,0) = (TT)1;
        return;
    }
    
    MatrixDecompTridiag<TT> tri(mavA);
    
    m_vdiag.copy_from(tri.diagonal());
    m_vsubdiag.copy_from(tri.sub_diagonal());
    
    Vec<TT> vgivens_c(size);
    Vec<TT> vgivens_s(size);
    
    if(compute_vectors)
        m_mvectors = tri.matrix_Q();
    
    core_priv_chop_small(m_vdiag.view(), m_vsubdiag.view());
    
    // reduce until diagonal

    int b = size-1;
    while(b>0)
    {
        if(m_vsubdiag(b-1)==(TT)0 || std::isnan(m_vsubdiag(b-1)))
        {
            b--;
            continue;
        }
        
        // find the largest unreduced block (a,b) starting from b and working backwards
        int a = b-1;
        while(a>0)
        {
            if(m_vsubdiag(a-1)==(TT)0)
                break;
            a--;
        }
        
        int nblock = b-a+1;
        
        VecView<TT> vdiagblock = m_vdiag.view().block(a, nblock);
        VecView<TT> vsubdiagblock = m_vsubdiag.view().block(a, nblock-1);
       
        // Apply QR reduction with implicit deflation to the unreduced block
        if(compute_vectors)
        {
            core_priv_qrstep(vdiagblock, vsubdiagblock, vgivens_c.view(), vgivens_s.view());
            
            // Apply givens rotation Gij(c,s) to Q
            for(int i=0; i<nblock-1; i++)
            {
                TT c = vgivens_c(i);
                TT s = vgivens_s(i);
                
                for(int k=0; k<size; k++)
                {
                    TT qki = m_mvectors(k,a+i);
                    TT qkj = m_mvectors(k,a+i+1);
                    
                    m_mvectors(k,a+i) = qki * c - qkj * s;
                    m_mvectors(k,a+i+1) = qki * s + qkj * c;
                }
            }
        }
        else
            core_priv_qrstep(vdiagblock, vsubdiagblock, vgivens_c.view(), vgivens_s.view());
        
        core_priv_chop_small(vsubdiagblock, vdiagblock);
    }
    
    // sometimes for (near) singular matrices eigenvalues can be tiny negatives due to float accuracy issues.
    for(int i=0; i<size; i++)
        if(m_vdiag(i)<(TT)0)
            m_vdiag(i) = (TT)0;
    
    // sort results
    for(int i=0; i<size-1; i++)
    {
        int k = i + core_block_location_of_max(m_vdiag.view().tail(size-i));
        std::swap(m_vdiag(i), m_vdiag(k));
        if(compute_vectors)
            core_block_exchange(m_mvectors.view().col(i), m_mvectors.view().col(k));
    }
    
    m_vectors = compute_vectors;
}

template <typename TT>
im::Mtx<TT> im::MatrixDecompEigenSymmetric<TT>::solve(MtxView<TT> const &mavy) const
{
    if(!m_vectors)
        IM_THROW_NO_INIT;
    
    IM_CHECK_VALID(mavy);
    int n = m_vdiag.rows();
    IM_CHECK_ARGS(mavy.rows()==n);
    
    Vec<TT> vtmp(n);
    Mtx<TT> mrtn(n,mavy.cols());
    
    for(int i=0; i<mavy.cols(); i++)
    {
        core_block_blas_gemv(vtmp.view(), m_mvectors.view().t(), mavy.col(i), (TT)1, (TT)0, TransMode_N);
        
        core_block_divide_elements(vtmp.view(), vtmp.view(), m_vdiag.view());
        
        core_block_blas_gemv(mrtn.view().col(i), m_mvectors.view(), vtmp.view(), (TT)1, (TT)0, TransMode_N);
    }
    
    return mrtn;
}

template <typename TT>
im::Vec<TT> im::MatrixDecompEigenSymmetric<TT>::solve(VecView<TT> const &vvy) const
{
    if(!m_vectors)
        IM_THROW_NO_INIT;
    
    IM_CHECK_VALID(vvy);
    int n = m_vdiag.rows();
    IM_CHECK_ARGS(vvy.rows()==n);
    
    Vec<TT> vtmp(n);
    Vec<TT> vrtn(n);
    
    core_block_blas_gemv(vtmp.view(), m_mvectors.view().t(), vvy, (TT)1, (TT)0, TransMode_N);
    
    core_block_divide_elements(vtmp.view(), vtmp.view(), m_vdiag.view());
    
    core_block_blas_gemv(vrtn.view(), m_mvectors.view(), vtmp.view(), (TT)1, (TT)0, TransMode_N);
    
    return vrtn;
}

// Computes V D^(1/2) V'
template <typename TT>
im::Mtx<TT> im::MatrixDecompEigenSymmetric<TT>::square_root() const
{
    if(!m_vectors)
        IM_THROW_NO_INIT;
    
    int n = m_vdiag.rows();
    Mtx<TT> mrtn(n, n);
    Mtx<TT> mtmp(n, n);
    Vec<TT> vtmp(n);
    
    core_block_sqrt(vtmp.view(), m_vdiag.view());
    
    core_block_multiply_rows(mtmp.view(), vtmp.view(), m_mvectors.view().t());
    
    core_block_blas_gemm(mrtn.view(), m_mvectors.view(), mtmp.view(), (TT)1, (TT)0, TransMode_N, TransMode_N);
    
    return mrtn;
}

// Computes V D^(-1/2) V'
template <typename TT>
im::Mtx<TT> im::MatrixDecompEigenSymmetric<TT>::inverse_square_root() const
{
    if(!m_vectors)
        IM_THROW_NO_INIT;
    
    int n = m_vdiag.rows();
    Mtx<TT> mrtn(n, n);
    Mtx<TT> mtmp(n, n);
    Vec<TT> vtmp(n);
    
    core_block_sqrt(vtmp.view(), m_vdiag.view());
    core_block_recip(vtmp.view(), vtmp.view());
    
    core_block_multiply_rows(mtmp.view(), vtmp.view(), m_mvectors.view().t());
    
    core_block_blas_gemm(mrtn.view(), m_mvectors.view(), mtmp.view(), (TT)1, (TT)0, TransMode_N, TransMode_N);
    
    return mrtn;
}

template class im::MatrixDecompEigenSymmetric<float>;
template class im::MatrixDecompEigenSymmetric<double>;

template <typename TT> void im::core_decomp_eigen_2x2(VecView<TT> vevals, MtxView<TT> mevecs, MtxView<TT> const &mA)
{
    IM_CHECK_VALID(vevals);
    IM_CHECK_VALID(mevecs);
    IM_CHECK_VALID(mA);
    IM_CHECK_VECTOR_SIZE(vevals, 2);
    IM_CHECK_MATRIX_SIZE(mA, 2, 2);
    IM_CHECK_MATRIX_SIZE(mevecs, 2, 2);

    TT a = mA(0,0);
    TT b = mA(0,1);
    TT c = mA(1,1);
    
    TT coefs[3];
    TT evec[2];

    coefs[0] = 1;
    coefs[1] = -a-c;
    coefs[2] = a*c - b*b;
    
    core_solve_quadratic_pos_real(vevals, VecView<TT>(3, 1, coefs));
    
    if(vevals(0)<=(TT)0)
    {
        mevecs = (TT)0;
        mevecs.diag() = (TT)1;
        return;
    }
    
    TT d0 = a - vevals(1);
    TT d1 = c - vevals(1);
    
    TT fb = std::abs(b);
    TT r0 = std::abs(d0) + fb;
    TT r1 = std::abs(d1) + fb;
    
    if(r1 > r0)
    {
        evec[0] = d1;
        evec[1] = -b;
    }
    else if(r0 > (TT)0)
    {
        evec[0] = b;
        evec[1] = -d0;
    }
    else
    {
        evec[0] = (TT)1;
        evec[1] = (TT)0;
    }
    
    TT sums = std::sqrt(evec[0]*evec[0]+evec[1]*evec[1]);
    evec[0]/=sums;
    evec[1]/=sums;
    
    mevecs(0,0) = -evec[1];
    mevecs(0,1) = evec[0];
    mevecs(1,0) = evec[0];
    mevecs(1,1) = evec[1];
    
}

template void im::core_decomp_eigen_2x2(VecView<float> vevals, MtxView<float> mevecs, MtxView<float> const &mA);
template void im::core_decomp_eigen_2x2(VecView<double> vevals, MtxView<double> mevecs, MtxView<double> const &mA);

