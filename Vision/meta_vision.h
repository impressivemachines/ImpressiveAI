//
//  meta_vision.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 3/7/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_meta_vision_h
#define Metaphor_meta_vision_h

namespace im
{
    // Color model used by images
    // Defines layout of multi-band images in memory
    enum ColorModel
    {
        ColorModelUndefined = 0,
        ColorModelGray,
        ColorModelRGBA,
        ColorModelCMYK,
    };
    
    // Integer point, size, and rect objects for images
    
    // point
    struct ImgPoint
    {
        ImgPoint() {}
        ImgPoint(int xx, int yy) : x(xx), y(yy) {}
        
        int x, y;
    };
    
    // size
    struct ImgSize
    {
        ImgSize() {}
        ImgSize(int ww, int hh) : width(ww), height(hh) {}
        
        bool is_empty() const { return width<1 || height<1; }
        ImgSize clip(ImgSize const &other) const { return ImgSize(std::min(width, other.width), std::min(height, other.height)); }
        void normalize() { width = abs(width); height = abs(height); }
        
        int width, height;
    };
    
    // rect
    struct ImgRect
    {
        ImgRect() {}
        ImgRect(int xx, int yy, int ww, int hh) : origin(ImgPoint(xx, yy)), size(ImgSize(ww, hh)) {}
        ImgRect(ImgPoint org, ImgSize sz) : origin(org), size(sz) {}
        
        bool is_empty() const { return size.is_empty(); }
        bool is_intersecting(ImgRect const &other) const;
        bool is_within(ImgRect const &other) const;
        ImgRect clip(ImgRect const &other) const;
        
        ImgPoint origin;
        ImgSize size;
    };
}

#include "Include/image.h"


#endif
