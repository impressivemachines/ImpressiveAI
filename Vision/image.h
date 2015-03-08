//
//  image.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/7/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_image_h
#define Metaphor_image_h

namespace im
{
    
    // Core information describing images of arbitrary packed raster format and type
    struct ImgInfo
    {
        ImgInfo()
        : size(ImgSize(0,0)), bands(0), elsize_bytes(0), stride_bytes(0), eltype(CoreTypeUndefined), color_model(ColorModelUndefined) {}
        
        ImgInfo(ImgSize sz, int b, size_t elsize, size_t stride, CoreType t, ColorModel cm = ColorModelUndefined)
        : size(sz), bands(b), elsize_bytes(elsize), stride_bytes(stride), eltype(t), color_model(cm) {}
        
        ImgSize size;
        int bands;
        size_t elsize_bytes; // size of one band
        size_t stride_bytes;
        CoreType eltype;
        ColorModel color_model;
    };
    
    // Non-templated generic image view for accessing images regardless of pixel type
    class GenericImgView
    {
    public:
        GenericImgView() : m_p(NULL), m_pixsize_bytes(0) {}
        GenericImgView(ImgInfo const &info, uint8_t const *p)
        : m_info(info), m_p(const_cast<uint8_t *>(p)), m_pixsize_bytes(info.bands * info.elsize_bytes) {}
        
        // wrap existing image data
        uint8_t *wrap(ImgInfo const &info, uint8_t const *p)
        {
            m_info = info;
            m_p = const_cast<uint8_t *>(p);
            m_pixsize_bytes = info.bands * info.elsize_bytes;
            
            // Stride has to be a multiple of elsize for things to work
            IM_CHECK_ARGS(info.stride_bytes == info.elsize_bytes * (info.stride_bytes / info.elsize_bytes));
            
            // pointer to next possible allocation location
            return m_p + stride_bytes() * height();
        }
        
        void reset()
        {
            m_info = ImgInfo();
            m_p = NULL;
            m_pixsize_bytes = 0;
        }
        
        // pointer to start of image
        uint8_t *byte_ptr() { return m_p; }
        uint8_t const *byte_ptr() const { return m_p; }
        
        // pointer to start of row
        uint8_t *byte_ptr(int y)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, m_info.size.height);
            return m_p + y * m_info.stride_bytes;
        }
        
        uint8_t const *byte_ptr(int y) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, m_info.size.height);
            return m_p + y * m_info.stride_bytes;
        }
        
        // pointer to start of pixel
        uint8_t *byte_ptr(int x, int y)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(x, 0, m_info.size.width);
            return m_p + y * m_info.stride_bytes + x * m_pixsize_bytes;
        }
        
        uint8_t const *byte_ptr(int x, int y) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(x, 0, m_info.size.width);
            return m_p + y * m_info.stride_bytes + x * m_pixsize_bytes;
        }
        
        // pointer to specific band
        uint8_t *byte_ptr(int x, int y, int b)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(b, 0, m_info.bands);
            return m_p + y * m_info.stride_bytes + x * m_pixsize_bytes + b * m_info.elsize_bytes;
        }
        
        uint8_t const *byte_ptr(int x, int y, int b) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(b, 0, m_info.bands);
            return m_p + y * m_info.stride_bytes + x * m_pixsize_bytes + b * m_info.elsize_bytes;
        }
        
        // return a sub block of the image
        GenericImgView block(int x, int y, int w, int h) const
        {
            IM_CHECK_ARGS(x>=0 && x+w <= m_info.size.width);
            IM_CHECK_ARGS(y>=0 && y+h <= m_info.size.height);
            
            ImgInfo info = m_info;
            info.size.width = w;
            info.size.width = h;
            return GenericImgView(info, byte_ptr(x,y));
        }
        
        GenericImgView block(ImgPoint pt, ImgSize sz) const
        {
            return block(pt.x, pt.y, sz.width, sz.height);
        }
        
        GenericImgView block(ImgRect rct) const
        {
            return block(rct.origin.x, rct.origin.y, rct.size.width, rct.size.height);
        }
        
        // clear the image to zero
        void clear();
        
        // copy from another view with exactly the same size and format
        void copy_from(GenericImgView const &gv);
        
        // paste from another view with any  source size and at position, including out of destination bounds
        void paste_from(GenericImgView const &gv, int destx, int desty);
        
        // get information about the image
        bool is_valid() const { return m_p!=NULL; }
        bool is_aligned_16() const { return (stride_bytes() & 0xf)==0 && IM_IS_ALIGNED_16(m_p); }
        ImgInfo info() const { return m_info; }
        int width() const { return m_info.size.width; }
        int height() const { return m_info.size.height; }
        size_t count() const { return (size_t)width() * (size_t)height(); }
        ImgSize size() const { return m_info.size; }
        ImgRect rect() const { return ImgRect(ImgPoint(0,0), m_info.size); }
        int bands() const { return m_info.bands; }
        size_t elsize_bytes() const { return m_info.elsize_bytes; }
        size_t pixsize_bytes() const { return m_pixsize_bytes; }
        size_t stride_bytes() const { return m_info.stride_bytes; }
        ColorModel color_model() const { return m_info.color_model; }
        CoreType type() const { return m_info.eltype; }
        
        // print the size of the image
        void print_size(bool cr = true, FILE *fp = stdout) const
        {
            IM_CHECK_NULL(fp);
            fprintf(fp, "(%d x %d)", width(), height());
            if(cr)
                fputc('\n', fp);
        }
        
    private:
        ImgInfo m_info; // format of the data
        size_t m_pixsize_bytes; // optimization
        uint8_t *m_p; // pointer to the data
    };
    
    // Templated image view that works with the specific element types supported by the library:
    
    // CoreTypeByte (uint8_t)
    // CoreTypeShort (int16_t)
    // CoreTypeInt (int)
    // CoreTypeFloat (float)
    // CoreTypeDouble (double)
    // CoreTypeComplexFloat (std::complex<float>)
    // CoreTypeComplexDouble (std::complex<double>)
    
    template <typename TT>
    class ImgView
    {
    public:
        ImgView() {}
        
        ImgView(GenericImgView const &gv) { wrap(gv); }
        
        ImgView(ImgInfo const &info, TT const *p) { wrap(info, p); }
        
        ImgView(int w, int h, int b, int row_stride, TT const *p, ColorModel cm = ColorModelUndefined)
            { wrap(w, h, b, row_stride, p, cm); }
        
        // wrap existing image data
        void wrap(GenericImgView const &gv) { wrap(gv.info(), (TT const *)gv.byte_ptr()); }
        
        // wrap existing image data
        TT *wrap(ImgInfo const &info, TT const *p)
        {
            IM_CHECK_ARGS(info.eltype==type());
            IM_CHECK_ARGS(info.elsize_bytes==sizeof(TT));
            IM_CHECK_ARGS(info.stride_bytes == info.elsize_bytes * (info.stride_bytes / info.elsize_bytes));
            
            IM_CHECK_LOWER_BOUNDS(info.size.width, 0);
            IM_CHECK_LOWER_BOUNDS(info.size.height, 0);
            IM_CHECK_LOWER_BOUNDS(info.bands, 0);
            
            m_size = info.size;
            m_bands = info.bands;
            m_stride = (int)(info.stride_bytes / sizeof(TT));
            m_color_model = info.color_model;
            m_p = const_cast<TT *>(p);
            
            return m_p + stride() * height();
        }
        
        // wrap existing image data
        // stride is in units of sizeof(TT) bytes
        // p can be NULL
        TT *wrap(int w, int h, int b, int row_stride, TT const *p, ColorModel cm = ColorModelUndefined)
        {
            IM_CHECK_LOWER_BOUNDS(w, 0);
            IM_CHECK_LOWER_BOUNDS(h, 0);
            IM_CHECK_LOWER_BOUNDS(b, 0);
            IM_CHECK_LOWER_BOUNDS(row_stride, 0);
            
            m_size.width = w;
            m_size.height = h;
            m_bands = b;
            m_stride = row_stride;
            m_color_model = cm;
            m_p = const_cast<TT *>(p);
            
            return m_p + stride() * height();
        }
        
        TT *ptr() { return m_p; }
        TT const *ptr() const { return m_p; }
        
        TT *ptr(int y)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, height());
            
            return m_p + y * m_stride;
        }
        
        TT const *ptr(int y) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, height());
            
            return m_p + y * m_stride;
        }
        
        TT *ptr(int x, int y)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, height());
            IM_DEBUG_ONLY_CHECK_BOUNDS(x, 0, width());
            
            return m_p + y * m_stride + x * m_bands;
        }
        
        TT const *ptr(int x, int y) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, height());
            IM_DEBUG_ONLY_CHECK_BOUNDS(x, 0, width());
            
            return m_p + y * m_stride + x * m_bands;
        }
        
        TT *ptr(int x, int y, int b)
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, height());
            IM_DEBUG_ONLY_CHECK_BOUNDS(x, 0, width());
            IM_DEBUG_ONLY_CHECK_BOUNDS(b, 0, bands());
            
            return m_p + y * m_stride + x * m_bands + b;
        }
        
        TT const *ptr(int x, int y, int b) const
        {
            IM_DEBUG_ONLY_CHECK_NULL(m_p);
            IM_DEBUG_ONLY_CHECK_BOUNDS(y, 0, height());
            IM_DEBUG_ONLY_CHECK_BOUNDS(x, 0, width());
            IM_DEBUG_ONLY_CHECK_BOUNDS(b, 0, bands());
            
            return m_p + y * m_stride + x * m_bands + b;
        }
        
        TT &at(int x, int y) { return *ptr(x,y); }
        TT const &at(int x, int y) const { return *ptr(x,y); }
        
        TT &at(int x, int y, int b) { return *ptr(x,y,b); }
        TT const &at(int x, int y, int b) const { return *ptr(x,y,b); }
        
        TT &operator()(int x, int y) { return *ptr(x,y); }
        TT const &operator()(int x, int y) const { return *ptr(x,y); }
        
        TT &operator()(int x, int y, int b) { return *ptr(x,y,b); }
        TT const &operator()(int x, int y, int b) const { return *ptr(x,y,b); }

        GenericImgView generic_view()
        {
            return GenericImgView(info(), (uint8_t *)m_p);
        }
        
        // matrix view for all bands
        MtxView<TT> matrix_view()
        {
            return MtxView<TT>(height(), width() * bands(), stride(), 1, m_p);
        }
        
        // matrix view for one band
        MtxView<TT> matrix_view(int band)
        {
            IM_CHECK_BOUNDS(band, 0, bands());
            return MtxView<TT>(height(), width(), stride(), bands(), m_p + band);
        }
        
        ImgView block(int x, int y, int w, int h) const
        {
            IM_CHECK_ARGS(x>=0 && x+w <= width());
            IM_CHECK_ARGS(y>=0 && y+h <= height());
            
            return ImgView(w, h, bands(), stride(), ptr(x,y), color_model());
        }
        
        ImgView block(ImgPoint pt, ImgSize sz) const { return block(pt.x, pt.y, sz.width, sz.height); }
        
        ImgView block(ImgRect rct) const { return block(rct.origin.x, rct.origin.y, rct.size.width, rct.size.height); }
        
        void clear() { fill((TT)0); }
        
        // fill all bands with a given value
        void fill(TT const &value)
        {
            IM_CHECK_ARGS(is_valid());
            matrix_view() = value;
        }
        
        ImgView &operator=(TT const &val)
        {
            fill(val);
            return *this;
        }
        
        // fill only the given band to value
        void fill(TT const &value, int band)
        {
            IM_CHECK_ARGS(is_valid());
            IM_CHECK_BOUNDS(band, 0, bands());
            matrix_view(band) = value;
        }
        
        // Fills all pixels with the multi-band data given
        void fill(TT const *ppixel);
        
        // copy from another view with exactly the same size and format
        void copy_from(ImgView const &iv);
        
        // paste from another view with any  source size and at position, including out of destination bounds
        void paste_from(ImgView const &iv, int destx, int desty);
        
        // get information about the image
        bool is_valid() const { return m_p!=NULL; }
        bool is_aligned_16() const { return (stride_bytes() & 0xf)==0 && IM_IS_ALIGNED_16(m_p); }
        ImgInfo info() const { return ImgInfo(m_size, m_bands, elsize_bytes(), stride_bytes(), type(), m_color_model); }
        int width() const { return m_size.width; }
        int height() const { return m_size.height; }
        size_t count() const { return (size_t)width() * (size_t)height(); }
        ImgSize size() const { return m_size; }
        ImgRect rect() const { return ImgRect(ImgPoint(0,0), m_size); }
        int bands() const { return m_bands; }
        size_t elsize_bytes() const { return sizeof(TT); }
        size_t pixsize_bytes() const { return m_bands * sizeof(TT); }
        size_t stride_bytes() const { return m_stride * sizeof(TT); }
        int stride() const { return m_stride; }
        ColorModel color_model() const { return m_color_model; }
        CoreType type() const { return (CoreType)TypeProperties<TT>::core_type_id; }
        
        // print the size of the image
        void print_size(bool cr = true, FILE *fp = stdout) const;
        
        // print out all the band values of a single pixel
        void print_pixel(int x, int y, bool cr = true, FILE *fp = stdout) const;
        
    private:
        ImgSize m_size;
        int m_bands;
        int m_stride;
        ColorModel m_color_model;
        TT *m_p;
    };

    // Templated image object that maintains shared pixel storage
    template <typename TT>
    class Img : public ImgView<TT>
    {
    public:
        Img() {}
        
    private:
    };
    
}


#endif