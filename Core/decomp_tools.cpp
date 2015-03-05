//
//  decomp_tools.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/31/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
TT im::core_compute_householder_vector(VecView<TT> vhh)
{
    IM_CHECK_VALID(vhh);
    
    int size = vhh.rows();
    if(size<2)
        return (TT)0;
    
    VecView<TT> vtail = vhh.tail(size-1);
    
    TT norm = std::sqrt(core_block_reduce_sum_squares(vtail));
    if(norm==(TT)0)
        return (TT)0;
    
    TT alpha = vhh(0);
    TT beta = core_hypot(alpha, norm);
    if(alpha >= (TT)0)
        beta = -beta;
    TT tau = (beta - alpha)/beta;
    TT s = alpha - beta;
    
    if(std::abs(s) > std::numeric_limits<TT>::min())
        core_block_scale(vtail, vtail, (TT)1 / s);
    else
    {
        // deal with denormals
        core_block_scale(vtail, vtail, std::numeric_limits<TT>::epsilon() / s);
        core_block_scale(vtail, vtail, (TT)1 / std::numeric_limits<TT>::epsilon());
    }
    
    vhh(0) = beta;
    
    return tau;
}

template float im::core_compute_householder_vector(VecView<float> v);
template double im::core_compute_householder_vector(VecView<double> v);

template <typename TT> void im::core_apply_householder_vector(MtxView<TT> mvA, TT const &tau, VecView<TT> const &vhh)
{
    IM_CHECK_VALID(mvA);
    IM_CHECK_VALID(vhh);
    
    // Applies a householder transform vector and tau to the given matrix
    // A = A - tau w v'
    
    if(tau==(TT)0)
        return;
    
    int rows = mvA.rows();
    int cols = mvA.cols();
    
    IM_CHECK_ARGS(rows>1);
    IM_CHECK_ARGS(rows==vhh.rows());
    
    int i;
    for(i=0; i<cols; i++)
    {
        VecView<TT> vAcol = mvA.col(i).tail(rows-1);
        VecView<TT> vHHsub = vhh.tail(rows-1);
        
        TT w = core_block_reduce_multiply_add(vAcol, vHHsub);
        w += mvA(0,i);
        TT scale = -tau * w;
        mvA(0,i) += scale;
        core_block_add_scaled(vAcol, vAcol, vHHsub, scale);
    }
}

template void im::core_apply_householder_vector(MtxView<float> mvA, float const &tau, VecView<float> const &vhh);
template void im::core_apply_householder_vector(MtxView<double> mvA, double const &tau, VecView<double> const &vhh);

template <typename TT> void im::core_apply_householder_identity(MtxView<TT> mvA, TT const &tau)
{
    IM_CHECK_VALID(mvA);
    
    if(tau==(TT)0)
    {
        mvA(0,0) = (TT)1;
        if(mvA.cols()>1)
            mvA.row(0).tail(mvA.cols()-1) = (TT)0;
        if(mvA.rows()>1)
            mvA.col(0).tail(mvA.rows()-1) = (TT)0;
        return;
    }
    
    if(mvA.rows()>1)
    {
        MtxView<TT> mvAsub = mvA.block(1, 0, mvA.rows()-1, mvA.cols());
        VecView<TT> vv0 = mvAsub.col(0);
        
        for(int j=1; j<mvA.cols(); j++)
        {
            VecView<TT> vvj = mvAsub.col(j);
            TT s = -tau * core_block_blas_dot(vv0, vvj);
            mvA(0,j) = s;
            core_block_blas_axpy(vvj, vv0, s);
        }
        
        core_block_scale(vv0, vv0, -tau);
    }
    
    mvA(0,0) = (TT)1 - tau;
}

template void im::core_apply_householder_identity(MtxView<float> mavA, float const &tau);
template void im::core_apply_householder_identity(MtxView<double> mavA, double const &tau);


