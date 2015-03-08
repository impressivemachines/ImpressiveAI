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
    size_t count = m_info.size.width * m_pixsize_bytes;
    for(int y=0; y<m_info.size.height; y++)
        memset(byte_ptr(y), 0, count);
}

void im::GenericImgView::copy_from(GenericImgView const &gv)
{
    IM_CHECK_ARGS(gv.width()==width());
    IM_CHECK_ARGS(gv.height()==height());
    IM_CHECK_ARGS(gv.bands()==bands());
    IM_CHECK_ARGS(gv.elsize_bytes()==elsize_bytes());
    IM_CHECK_ARGS(gv.type()==type());
    
    size_t count = width() * pixsize_bytes();
    for(int y=0; y<m_info.size.height; y++)
        memcpy(byte_ptr(y), gv.byte_ptr(y), count);
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

template <typename TT>
void im::ImgView<TT>::fill(TT const *ppixel)
{
    IM_CHECK_ARGS(is_valid());
    IM_CHECK_NULL(ppixel);
    
    for(int x=0; x<width(); x++)
        memcpy(ptr(x,0), ppixel, pixsize_bytes());
    for(int y=1; y<height(); y++)
        memcpy(ptr(y), ptr(0), pixsize_bytes() * width());
}

template <typename TT>
void im::ImgView<TT>::copy_from(ImgView const &iv)
{
    IM_CHECK_ARGS(iv.width()==width());
    IM_CHECK_ARGS(iv.height()==height());
    IM_CHECK_ARGS(iv.bands()==bands());
    
    size_t count = width() * pixsize_bytes();
    for(int y=0; y<m_size.height; y++)
        memcpy(ptr(y), iv.ptr(y), count);
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

