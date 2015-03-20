//
//  image.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/7/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"
#include "meta_vision.h"


void im::GenericImgView::clear()
{
    size_t row_bytes = width() * pixsize_bytes();
    for(int y=0; y<m_info.size.height; y++)
        memset(ptr(y), 0, row_bytes);
}

void im::GenericImgView::copy_from(GenericImgView const &gv)
{
    IM_CHECK_ARGS(gv.width()==width());
    IM_CHECK_ARGS(gv.height()==height());
    IM_CHECK_ARGS(gv.bands()==bands());
    IM_CHECK_ARGS(gv.elsize_bytes()==elsize_bytes());
    IM_CHECK_ARGS(gv.type()==type());
    
    size_t row_bytes = width() * pixsize_bytes();
    for(int y=0; y<m_info.size.height; y++)
        memcpy(ptr(y), gv.ptr(y), row_bytes);
}

// Source and destination can be any size and position, including out of destination bounds
void im::GenericImgView::paste_from(GenericImgView const &gv, int destx, int desty)
{
    IM_CHECK_ARGS(gv.bands()==bands());
    IM_CHECK_ARGS(gv.elsize_bytes()==elsize_bytes());
    IM_CHECK_ARGS(gv.type()==type());
    
    ImgRect rctdest = ImgRect(destx, desty, gv.width(), gv.height());
    ImgRect rctcopydest = rctdest.clip(rect());
    if(!rctcopydest.is_empty())
    {
        ImgRect rctcopysrc = rctcopydest;
        rctcopysrc.origin.x -= destx;
        rctcopysrc.origin.y -= desty;
        block(rctcopydest).copy_from(gv.block(rctcopysrc));
    }
}
/*

template <typename TT>
void im::ImgView<TT>::fill(TT const *ppixel)
{
    IM_CHECK_ARGS(is_valid());
    IM_CHECK_NULL(ppixel);
    
    for(int x=0; x<width(); x++)
        memcpy(ptr(x,0), ppixel, pixsize_bytes());
    
    size_t row_bytes = width() * pixsize_bytes();
    for(int y=1; y<height(); y++)
        memcpy(ptr(y), ptr(0), row_bytes);
}

template <typename TT>
void im::ImgView<TT>::copy_from(ImgView const &iv)
{
    IM_CHECK_ARGS(iv.width()==width());
    IM_CHECK_ARGS(iv.height()==height());
    IM_CHECK_ARGS(iv.bands()==bands());
    
    size_t row_bytes = width() * pixsize_bytes();
    for(int y=0; y<m_size.height; y++)
        memcpy(ptr(y), iv.ptr(y), row_bytes);
}

// Source and destination can be any size and position, including out of destination bounds
template <typename TT>
void im::ImgView<TT>::paste_from(ImgView const &iv, int destx, int desty)
{
    IM_CHECK_ARGS(iv.bands()==bands());
    
    ImgRect rctdest = ImgRect(destx, desty, iv.width(), iv.height());
    ImgRect rctcopydest = rctdest.clip(rect());
    if(!rctcopydest.is_empty())
    {
        ImgRect rctcopysrc = rctcopydest;
        rctcopysrc.origin.x -= destx;
        rctcopysrc.origin.y -= desty;
        block(rctcopydest).copy_from(iv.block(rctcopysrc));
    }
}

template <typename TT>
void im::ImgView<TT>::print_size(bool cr, FILE *fp) const
{
    IM_CHECK_NULL(fp);
    fprintf(fp, "(%d x %d)", width(), height());
    if(cr)
        fputc('\n', fp);
}

template <typename TT>
void im::ImgView<TT>::print_pixel(int x, int y, bool cr, FILE *fp) const
{
    IM_CHECK_NULL(fp);
    fprintf(fp, "(%d,%d)=[", x, y);
    for(int b=0; b<bands(); b++)
    {
        core_print_value(fp, at(x,y,b));
        if(b+1 != bands())
            fputc(' ', fp);
    }
    fputc(']', fp);
    if(cr)
        fputc('\n', fp);
}

// Instantiate classes
template class im::ImgView<uint8_t>;
template class im::ImgView<int16_t>;
template class im::ImgView<int>;
template class im::ImgView<float>;
template class im::ImgView<double>;
template class im::ImgView< std::complex<float> >;
template class im::ImgView< std::complex<double> >;

template <typename TT>
void im::Img<TT>::resize(int w, int h, int b, ColorModel cm)
{
    if(b<1)
        b = std::max(bands(),1);
    
    IM_CHECK_LOWER_BOUNDS(w, 0);
    IM_CHECK_LOWER_BOUNDS(h, 0);
    
    if(!(m_mem && m_mem.unique()))
        m_mem.reset(new MemoryBlock);
    
    // resize
    size_t rowbytes = w * b * sizeof(TT);
    
    // pad row to multiple of SSE size
    int pad = 0;
    while(!IM_IS_SIZED_SSE(rowbytes))
    {
        rowbytes += sizeof(TT);
        pad++;
        IM_CHECK(pad<IM_SSE_ALIGN);
    }
    
    int rowstride = (int)(rowbytes / sizeof(TT));
    
    m_mem->resize(rowstride * h, IM_SSE_ALIGN);
    
    m_view.wrap(w, h, b, rowstride, (TT *)m_mem->ptr(), cm);
}

template <typename TT>
void im::Img<TT>::stop_sharing()
{
    if(!m_mem || !m_mem.unique())
    {
        // Ensure simple array
        MemoryBlock *p = new MemoryBlock;
        
        // resize
        size_t rowbytes = width() * bands() * sizeof(TT);
        
        // pad row to multiple of SSE size
        int pad = 0;
        while(!IM_IS_SIZED_SSE(rowbytes))
        {
            rowbytes += sizeof(TT);
            pad++;
            IM_CHECK(pad<IM_SSE_ALIGN);
        }
        
        int rowstride = (int)(rowbytes / sizeof(TT));
        
        p->resize(rowstride * height(), IM_SSE_ALIGN);
        
        // copy over the data from the current view
        if(m_view.is_valid())
            ImgView<TT>(info(), (TT *)p->ptr()).copy_from(m_view);
        
        // switch to the new memory store
        m_mem.reset(p);
        m_view.wrap(info(), (TT *)p->ptr());
    }
}

// Instantiate classes
template class im::Img<uint8_t>;
template class im::Img<int16_t>;
template class im::Img<int>;
template class im::Img<float>;
template class im::Img<double>;
template class im::Img< std::complex<float> >;
template class im::Img< std::complex<double> >;
*/
