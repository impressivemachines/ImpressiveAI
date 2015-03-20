//
//  block_sort.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/25/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

template <typename TT> int core_priv_qsort_compare_ascending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( *(TT const *)p1 );
    double v2 = (double)im::core_real( *(TT const *)p2 );
    
    return v1 > v2 ? 1 :(v1 < v2 ? -1 : 0);
}

template <typename TT> int core_priv_qsort_compare_descending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( *(TT const *)p1 );
    double v2 = (double)im::core_real( *(TT const *)p2 );
    
    return v1 > v2 ? -1 :(v1 < v2 ? 1 : 0);
}

template <typename TT> int core_priv_qsort_compare_abs_ascending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( std::abs(*(TT const *)p1) );
    double v2 = (double)im::core_real( std::abs(*(TT const *)p2) );
    
    return v1 > v2 ? 1 :(v1 < v2 ? -1 : 0);
}

template <typename TT> int core_priv_qsort_compare_abs_descending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( std::abs(*(TT const *)p1) );
    double v2 = (double)im::core_real( std::abs(*(TT const *)p2) );
    
    return v1 > v2 ? -1 :(v1 < v2 ? 1 : 0);
}

template <typename TT> void im::core_vector_sort(std::vector<TT> &v, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        qsort((void *)v.data(), v.size(), sizeof(TT), core_priv_qsort_compare_ascending<TT>);
    else
        qsort((void *)v.data(), v.size(), sizeof(TT), core_priv_qsort_compare_descending<TT>);
}

#define INST(TT) template void im::core_vector_sort(std::vector<TT> &v, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST


template <typename TT> void im::core_vector_sort_abs(std::vector<TT> &v, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        qsort((void *)v.data(), v.size(), sizeof(TT), core_priv_qsort_compare_abs_ascending<TT>);
    else
        qsort((void *)v.data(), v.size(), sizeof(TT), core_priv_qsort_compare_abs_descending<TT>);
}

#define INST(TT) template void im::core_vector_sort_abs(std::vector<TT> &v, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> int core_priv_qsort_pointer_compare_ascending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( **(TT const **)p1 );
    double v2 = (double)im::core_real( **(TT const **)p2 );
    
    return v1 > v2 ? 1 :(v1 < v2 ? -1 : 0);
}

template <typename TT> int core_priv_qsort_pointer_compare_descending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( **(TT const **)p1 );
    double v2 = (double)im::core_real( **(TT const **)p2 );
    
    return v1 > v2 ? -1 :(v1 < v2 ? 1 : 0);
}

template <typename TT> int core_priv_qsort_pointer_compare_abs_ascending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( std::abs(**(TT const **)p1) );
    double v2 = (double)im::core_real( std::abs(**(TT const **)p2) );
    
    return v1 > v2 ? 1 :(v1 < v2 ? -1 : 0);
}

template <typename TT> int core_priv_qsort_pointer_compare_abs_descending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( std::abs(**(TT const **)p1) );
    double v2 = (double)im::core_real( std::abs(**(TT const **)p2) );
    
    return v1 > v2 ? -1 :(v1 < v2 ? 1 : 0);
}

//

template <typename TT> void core_priv_block_sort(im::VecView<TT> vdst, int (*compare_func)(void const *p1, void const *p2))
{
    IM_CHECK_VALID(vdst);
    
    int count = vdst.count();
    
    std::vector<TT const *> pointer_list;
    pointer_list.resize(count);
    
    // It is necessary to sort a pointer list because the row strides can be anything, and qsort expects sequential elements
    TT const **pplist = pointer_list.data();
    
    for(int row = 0; row<vdst.rows(); row++)
            *pplist++ = vdst.ptr(row);
    
    qsort((void *)pointer_list.data(), count, sizeof(TT const *), compare_func);
    
    std::vector<TT> vsorted;
    vsorted.resize(count);
    
    TT *psorted = vsorted.data();
    
    for(int i=0; i<count; i++)
        *psorted++ = *(pointer_list[i]);
    
    vdst.copy_from(im::VecView<TT>(vsorted));
}

template <typename TT> void im::core_block_sort(VecView<TT> vdst, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(vdst, core_priv_qsort_pointer_compare_ascending<TT>);
    else
        core_priv_block_sort(vdst, core_priv_qsort_pointer_compare_descending<TT>);
}

#define INST(TT) template void im::core_block_sort(VecView<TT> vdst, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sort_abs(VecView<TT> vdst, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(vdst, core_priv_qsort_pointer_compare_abs_ascending<TT>);
    else
        core_priv_block_sort(vdst, core_priv_qsort_pointer_compare_abs_descending<TT>);
}

#define INST(TT) template void im::core_block_sort_abs(VecView<TT> vdst, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void core_priv_block_sort(im::MtxView<TT> mdst, int (*compare_func)(void const *p1, void const *p2))
{
    IM_CHECK_VALID(mdst);
    
    int count = mdst.count();
    
    std::vector<TT const *> pointer_list;
    pointer_list.resize(count);
    
    // It is necessary to sort a pointer list because the row and column strides can be anything, and qsort expects sequential elements
    TT const **pplist = pointer_list.data();
    
    for(int row = 0; row<mdst.rows(); row++)
        for(int col = 0; col<mdst.cols(); col++)
            *pplist++ = mdst.ptr(row, col);
    
    qsort((void *)pointer_list.data(), count, sizeof(TT const *), compare_func);
    
    std::vector<TT> vsorted;
    vsorted.resize(count);
    
    TT *psorted = vsorted.data();
    
    for(int i=0; i<count; i++)
        *psorted++ = *(pointer_list[i]);
    
    mdst.copy_from(im::VecView<TT>(vsorted).rm_2d_matrix_view(mdst.cols()));
}

template <typename TT> void im::core_block_sort(MtxView<TT> mdst, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(mdst, core_priv_qsort_pointer_compare_ascending<TT>);
    else
        core_priv_block_sort(mdst, core_priv_qsort_pointer_compare_descending<TT>);
}

#define INST(TT) template void im::core_block_sort(MtxView<TT> mdst, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sort_abs(MtxView<TT> mdst, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(mdst, core_priv_qsort_pointer_compare_abs_ascending<TT>);
    else
        core_priv_block_sort(mdst, core_priv_qsort_pointer_compare_abs_descending<TT>);
}

#define INST(TT) template void im::core_block_sort_abs(MtxView<TT> mdst, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT>
struct QsortIndexTuple
{
    im::MtxLoc loc;
    TT const *pdata;
};

template <typename TT> int core_priv_qsort_tuple_compare_ascending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real(*((QsortIndexTuple<TT> const *)p1)->pdata);
    double v2 = (double)im::core_real(*((QsortIndexTuple<TT> const *)p2)->pdata);
    
    return v1 > v2 ? 1 :(v1 < v2 ? -1 : 0);
}

template <typename TT> int core_priv_qsort_tuple_compare_descending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real(*((QsortIndexTuple<TT> const *)p1)->pdata);
    double v2 = (double)im::core_real(*((QsortIndexTuple<TT> const *)p2)->pdata);
    
    return v1 > v2 ? -1 :(v1 < v2 ? 1 : 0);
}

template <typename TT> int core_priv_qsort_tuple_compare_abs_ascending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( std::abs(*((QsortIndexTuple<TT> const *)p1)->pdata) );
    double v2 = (double)im::core_real( std::abs(*((QsortIndexTuple<TT> const *)p2)->pdata) );
    
    return v1 > v2 ? 1 :(v1 < v2 ? -1 : 0);
}

template <typename TT> int core_priv_qsort_tuple_compare_abs_descending(void const *p1, void const *p2)
{
    double v1 = (double)im::core_real( std::abs(*((QsortIndexTuple<TT> const *)p1)->pdata) );
    double v2 = (double)im::core_real( std::abs(*((QsortIndexTuple<TT> const *)p2)->pdata) );
    
    return v1 > v2 ? -1 :(v1 < v2 ? 1 : 0);
}

//

template <typename TT> void core_priv_block_sort(std::vector<int> &index_list, im::VecView<TT> const &vsrc, int (*compare_func)(void const *p1, void const *p2))
{
    IM_CHECK_VALID(vsrc);
    
    int count = vsrc.count();
    
    index_list.resize(count);
    
    std::vector<QsortIndexTuple<TT> > tuple_list;
    tuple_list.resize(count);
    
    for(int row = 0; row<vsrc.rows(); row++)
    {
        QsortIndexTuple<TT> &tuple = tuple_list[row];
        tuple.loc.row = row;
        tuple.loc.col = 0;
        tuple.pdata = vsrc.ptr(row);
    }
    
    qsort((void *)tuple_list.data(), count, sizeof(QsortIndexTuple<TT>), compare_func);
    
    for(int i=0; i<count; i++)
        index_list[i] = tuple_list[i].loc.row;
}

template <typename TT> void im::core_block_sort(std::vector<int> &index_list, VecView<TT> const &vsrc, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(index_list, vsrc, core_priv_qsort_tuple_compare_ascending<TT>);
    else
        core_priv_block_sort(index_list, vsrc, core_priv_qsort_tuple_compare_descending<TT>);
}

#define INST(TT) template void im::core_block_sort(std::vector<int> &index_list, VecView<TT> const &vsrc, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sort_abs(std::vector<int> &index_list, VecView<TT> const &vsrc, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(index_list, vsrc, core_priv_qsort_tuple_compare_abs_ascending<TT>);
    else
        core_priv_block_sort(index_list, vsrc, core_priv_qsort_tuple_compare_abs_descending<TT>);
}

#define INST(TT) template void im::core_block_sort_abs(std::vector<int> &index_list, VecView<TT> const &vsrc, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> void core_priv_block_sort(std::vector<im::MtxLoc> &index_list, im::MtxView<TT> const &msrc, int (*compare_func)(void const *p1, void const *p2))
{
    IM_CHECK_VALID(msrc);
    
    int count = msrc.count();
    
    index_list.resize(count);
    
    std::vector<QsortIndexTuple<TT> > tuple_list;
    tuple_list.resize(count);
    
    int cols = msrc.cols();
    for(int row = 0; row<msrc.rows(); row++)
        for(int col = 0; col<cols; col++)
        {
            QsortIndexTuple<TT> &tuple = tuple_list[row * cols + col];
            tuple.loc.row = row;
            tuple.loc.col = col;
            tuple.pdata = msrc.ptr(row,col);
        }
    
    qsort((void *)tuple_list.data(), count, sizeof(QsortIndexTuple<TT>), compare_func);
    
    for(int i=0; i<count; i++)
        index_list[i] = tuple_list[i].loc;
}

template <typename TT> void im::core_block_sort(std::vector<MtxLoc> &index_list, MtxView<TT> const &msrc, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(index_list, msrc, core_priv_qsort_tuple_compare_ascending<TT>);
    else
        core_priv_block_sort(index_list, msrc, core_priv_qsort_tuple_compare_descending<TT>);
}

#define INST(TT) template void im::core_block_sort(std::vector<MtxLoc> &index_list, MtxView<TT> const &msrc, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_block_sort_abs(std::vector<MtxLoc> &index_list, MtxView<TT> const &msrc, SortDirection direction)
{
    if(direction==im::SortDirectionAscending)
        core_priv_block_sort(index_list, msrc, core_priv_qsort_tuple_compare_abs_ascending<TT>);
    else
        core_priv_block_sort(index_list, msrc, core_priv_qsort_tuple_compare_abs_descending<TT>);
}

#define INST(TT) template void im::core_block_sort_abs(std::vector<MtxLoc> &index_list, MtxView<TT> const &msrc, SortDirection direction)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> int im::core_block_location_of_max(VecView<TT> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    
    int rows = vsrc.rows();
    int best = 0;
    
    if(rows==0)
        return best;
    
    TT bestval = vsrc(0);
    int rowstride = vsrc.row_stride();
    TT const *p = vsrc.ptr();
    
    for(int i=0; i<rows; i++, p += rowstride)
    {
        TT val = *p;
        if(core_real(val) > core_real(bestval))
        {
            bestval = val;
            best = i;
        }
    }
    
    return best;
}

#define INST(TT) template int im::core_block_location_of_max(VecView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> int im::core_block_location_of_max_abs(VecView<TT> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    
    int rows = vsrc.rows();
    int best = 0;
    
    if(rows==0)
        return best;
    
    TT bestval = std::abs(vsrc(0));
    int rowstride = vsrc.row_stride();
    TT const *p = vsrc.ptr();
    
    for(int i=0; i<rows; i++, p += rowstride)
    {
        TT val = std::abs(*p);
        if(core_real(val) > core_real(bestval))
        {
            bestval = val;
            best = i;
        }
    }
    
    return best;
}

#define INST(TT) template int im::core_block_location_of_max_abs(VecView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> int im::core_block_location_of_min(VecView<TT> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    
    int rows = vsrc.rows();
    int best = 0;
    
    if(rows==0)
        return best;
    
    TT bestval = vsrc(0);
    int rowstride = vsrc.row_stride();
    TT const *p = vsrc.ptr();
    
    for(int i=0; i<rows; i++, p += rowstride)
    {
        TT val = *p;
        if(core_real(val) < core_real(bestval))
        {
            bestval = val;
            best = i;
        }
    }
    
    return best;
}

#define INST(TT) template int im::core_block_location_of_min(VecView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> int im::core_block_location_of_min_abs(VecView<TT> const &vsrc)
{
    IM_CHECK_VALID(vsrc);
    
    int rows = vsrc.rows();
    int best = 0;
    
    if(rows==0)
        return best;
    
    TT bestval = std::abs(vsrc(0));
    int rowstride = vsrc.row_stride();
    TT const *p = vsrc.ptr();
    
    for(int i=0; i<rows; i++, p += rowstride)
    {
        TT val = std::abs(*p);
        if(core_real(val) < core_real(bestval))
        {
            bestval = val;
            best = i;
        }
    }
    
    return best;
}

#define INST(TT) template int im::core_block_location_of_min_abs(VecView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

template <typename TT> im::MtxLoc im::core_block_location_of_max(MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(msrc);
    
    int rows = msrc.rows();
    int cols = msrc.cols();
    
    MtxLoc best(0,0);
    
    if(rows==0 || cols==0)
        return best;
    
    TT bestval = msrc(0,0);
    
    int colstride = msrc.col_stride();
    
    for(int i=0; i<rows; i++)
    {
        TT const *p = msrc.ptr(i,0);
        for(int j=0; j<cols; j++, p += colstride)
        {
            TT val = *p;
            if(core_real(val) > core_real(bestval))
            {
                bestval = val;
                best.row = i;
                best.col = j;
            }
        }
    }
    
    return best;
}

#define INST(TT) template im::MtxLoc im::core_block_location_of_max(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> im::MtxLoc im::core_block_location_of_max_abs(MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(msrc);
    
    int rows = msrc.rows();
    int cols = msrc.cols();
    
    MtxLoc best(0,0);
    
    if(rows==0 || cols==0)
        return best;
    
    TT bestval = std::abs(msrc(0,0));
    
    int colstride = msrc.col_stride();
    
    for(int i=0; i<rows; i++)
    {
        TT const *p = msrc.ptr(i,0);
        for(int j=0; j<cols; j++, p += colstride)
        {
            TT val = std::abs(*p);
            if(core_real(val) > core_real(bestval))
            {
                bestval = val;
                best.row = i;
                best.col = j;
            }
        }
    }
    
    return best;
}

#define INST(TT) template im::MtxLoc im::core_block_location_of_max_abs(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> im::MtxLoc im::core_block_location_of_min(MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(msrc);
    
    int rows = msrc.rows();
    int cols = msrc.cols();
    
    MtxLoc best(0,0);
    
    if(rows==0 || cols==0)
        return best;
    
    TT bestval = msrc(0,0);
    
    int colstride = msrc.col_stride();
    
    for(int i=0; i<rows; i++)
    {
        TT const *p = msrc.ptr(i,0);
        for(int j=0; j<cols; j++, p += colstride)
        {
            TT val = *p;
            if(core_real(val) < core_real(bestval))
            {
                bestval = val;
                best.row = i;
                best.col = j;
            }
        }
    }
    
    return best;
}

#define INST(TT) template im::MtxLoc im::core_block_location_of_min(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> im::MtxLoc im::core_block_location_of_min_abs(MtxView<TT> const &msrc)
{
    IM_CHECK_VALID(msrc);
    
    int rows = msrc.rows();
    int cols = msrc.cols();
    
    MtxLoc best(0,0);
    
    if(rows==0 || cols==0)
        return best;
    
    TT bestval = std::abs(msrc(0,0));
    
    int colstride = msrc.col_stride();
    
    for(int i=0; i<rows; i++)
    {
        TT const *p = msrc.ptr(i,0);
        for(int j=0; j<cols; j++, p += colstride)
        {
            TT val = std::abs(*p);
            if(core_real(val) < core_real(bestval))
            {
                bestval = val;
                best.row = i;
                best.col = j;
            }
        }
    }
    
    return best;
}

#define INST(TT) template im::MtxLoc im::core_block_location_of_min_abs(MtxView<TT> const &msrc)
INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST
