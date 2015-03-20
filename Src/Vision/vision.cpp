//
//  vision.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/7/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"
#include "meta_vision.h"

bool im::ImgRect::is_intersecting(ImgRect const &other) const
{
    int x0 = std::max(origin.x, other.origin.x);
    int x1 = std::min(origin.x + size.width, other.origin.x + other.size.width);
    if(x1 <= x0)
        return false;
    
    int y0 = std::max(origin.y, other.origin.y);
    int y1 = std::min(origin.y + size.height, other.origin.y + other.size.height);
    if(y1 <= y0)
        return false;
    
    return true;
}

bool im::ImgRect::is_within(ImgRect const &other) const
{
    if(origin.x < other.origin.x || origin.x + size.width > other.origin.x + other.size.width)
        return false;
    if(origin.y < other.origin.y || origin.y + size.height > other.origin.y + other.size.height)
        return false;
    
    return true;
}

im::ImgRect im::ImgRect::clip(ImgRect const &other) const
{
    int x0 = std::max(origin.x, other.origin.x);
    int x1 = std::min(origin.x + size.width, other.origin.x + other.size.width);
    int w = x1 - x0;
    int y0 = std::max(origin.y, other.origin.y);
    int y1 = std::min(origin.y + size.height, other.origin.y + other.size.height);
    int h = y1 - y0;
    
    if(w<1 || h<1)
        return ImgRect(origin, ImgSize(0,0));
    else
        return ImgRect(x0, y0, w, h);
}

