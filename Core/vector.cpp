//
//  vector.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/19/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

// Resizes the vec and loses any data
template <typename TT>
void im::Vec<TT>::resize(int nrows)
{
    IM_CHECK_LOWER_BOUNDS(nrows, 0);
    
    if(m_mem && m_mem.unique())
        m_mem->resize(nrows);
    else
        m_mem.reset(new std::vector<TT>(nrows));
    
    m_view.wrap(nrows, 1, m_mem->data());
}

// Stops sharing the vec by creating own private copy and guarentees packed row major layout
template <typename TT>
void im::Vec<TT>::stop_sharing()
{
    if(!m_mem || !m_mem.unique() || row_stride()!=1)
    {
        // Ensure simple array
        std::vector<TT> *pv = new std::vector<TT>(count());
        if(m_view.is_valid())
            core_block_copy(VecView<TT>(rows(), 1, pv->data()), m_view);
        m_mem.reset(pv);
        m_view.wrap(rows(), 1, m_mem->data());
    }
}

// Add an new row at the bottom of the vec and sets its element(s) to the given value
template <typename TT>
void im::Vec<TT>::push_back(TT const &f)
{
    stop_sharing(); // in case vec does not own memory
    
    if(rows()==0)
    {
        // Empty vec
        resize(1);
        at(0) = f;
        return;
    }

    int last = rows();
    m_mem->resize(last + 1);
    m_view.wrap(last + 1, 1, m_mem->data());
    at(last) = f;
}

// Remove a row from the bottom of the matrix
template <typename TT>
void im::Vec<TT>::pop_back()
{
    if(rows()>0)
    {
        stop_sharing(); // in case mtx does not own memory
        
        m_mem->resize(rows()-1);
        m_view.wrap(rows()-1, 1, m_mem->data());
    }
}

template <typename TT>
void im::Vec<TT>::random_uniform(Rand &rnd, TT const &low, TT const &high)
{
    for(int row=0; row<rows(); row++)
        at(row) = (TT)rnd.uniform_real(core_real(low), core_real(high));
}

template <typename TT>
void im::Vec<TT>::random_gaussian(Rand &rnd, TT const &mean, TT const &stddev)
{
    for(int row=0; row<rows(); row++)
        at(row) = stddev * (TT)rnd.gauss() + mean;
}

template <typename TT>
TT im::Vec<TT>::sample_bicubic(float row) const
{
    int nrows = rows();
    
    if(row<0 || row>nrows-1)
        return (TT)0;
    
    int r1 = (int)std::floor(row);
    int r0 = r1 - 1;
    int r2 = r1 + 1;
    int r3 = r1 + 2;
    
    if(r1<1 || r1>nrows-3)
    {
        r0 = core_range_limit(r0, 0, nrows-1);
        r2 = core_range_limit(r2, 0, nrows-1);
        r3 = core_range_limit(r3, 0, nrows-1);
    }
    
    float ar = row - r1;
    float mar = (1-ar);
    
    TT krow0 = (TT)((1/6.0f) * (mar * mar - 1.0f) * mar);
    TT krow1 = (TT)(0.5f * (ar * mar + 2.0f) * mar);
    TT krow2 = (TT)(0.5f * (mar * ar + 2.0f) * ar);
    TT krow3 = (TT)((1/6.0f) * (ar * ar - 1.0f) * ar);
    
    return krow0 * at(r0) + krow1 * at(r1) + krow2 * at(r2) + krow3 * at(r3);
}

template <typename TT>
TT im::Vec<TT>::sample_bilinear(float row) const
{
    if(row<0 || row>rows()-1)
        return (TT)0;
    
    int r0 = (int)std::floor(row);
    int r1 = std::min(r0+1, rows()-1);
    
    float ar = row - r0;
    
    TT v00 = at(r0);
    TT v01 = at(r1);
    
    return v00 + TT(ar) * (v01 - v00);
}

// Returns index list for values matching comparison
template <typename TT>
std::vector<int> im::Vec<TT>::select_equal_to(TT const &val)
{
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(at(i)==val)
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_greater_than(TT const &val)
{
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(core_real(at(i))>core_real(val))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_less_than(TT const &val)
{
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(core_real(at(i))<core_real(val))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_greater_than_abs(TT const &val)
{
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(std::abs(at(i)) > std::abs(val))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_less_than_abs(TT const &val)
{
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(std::abs(at(i)) < std::abs(val))
            index.push_back(i);
    
    return index;
}

// Returns index list for values matching comparison against elements of equal size vec

template <typename TT>
std::vector<int> im::Vec<TT>::select_equal_to(Vec const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_VECTOR_SIZES_MATCH(m_view, src);
    
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(at(i)==src(i))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_greater_than(Vec const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_VECTOR_SIZES_MATCH(m_view, src);
    
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(core_real(at(i)) > core_real(src(i)))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_less_than(Vec const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_VECTOR_SIZES_MATCH(m_view, src);
    
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(core_real(at(i)) < core_real(src(i)))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_greater_than_abs(Vec const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_VECTOR_SIZES_MATCH(m_view, src);
    
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(std::abs(at(i)) > std::abs(src(i)))
            index.push_back(i);
    
    return index;
}

template <typename TT>
std::vector<int> im::Vec<TT>::select_less_than_abs(Vec const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_VECTOR_SIZES_MATCH(m_view, src);
    
    std::vector<int> index;
    for(int i=0; i<rows(); i++)
        if(std::abs(at(i)) < std::abs(src(i)))
            index.push_back(i);
    
    return index;
}

// Instantiate classes
template class im::Vec<float>;
template class im::Vec<double>;
template class im::Vec<im::Cf>;
template class im::Vec<im::Cd>;
