//
//  utils.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"


template <typename TT> void im::core_sort_2(TT &a, TT &b, SortDirection dir)
{
    if(dir==im::SortDirectionAscending)
    {
        if(b < a)
            std::swap(a, b);
    }
    else
    {
        if(b > a)
            std::swap(a, b);
    }
}

template void im::core_sort_2(int &a, int &b, SortDirection dir);
template void im::core_sort_2(float &a, float &b, SortDirection dir);
template void im::core_sort_2(double &a, double &b, SortDirection dir);

template <typename TT> void im::core_sort_3(TT &a, TT &b, TT &c, SortDirection dir)
{
    if(dir==im::SortDirectionAscending)
    {
        if(a > c)
            std::swap(a, c);
        
        if(a > b)
            std::swap(a, b);
        
        if(b > c)
            std::swap(b, c);
    }
    else
    {
        if(a < c)
            std::swap(a, c);
        
        if(a < b)
            std::swap(a, b);
        
        if(b < c)
            std::swap(b, c);
    }
}

template void im::core_sort_3(int &a, int &b, int &c, SortDirection dir);
template void im::core_sort_3(float &a, float &b, float &c, SortDirection dir);
template void im::core_sort_3(double &a, double &b, double &c, SortDirection dir);

template <typename TT> void im::core_sort_4(TT &a, TT &b, TT &c, TT &d, SortDirection dir)
{
    if(dir==im::SortDirectionAscending)
    {
        if(a > b)
            std::swap(a, b);
        
        if(c > d)
            std::swap(c, d);
        
        if(a > c)
            std::swap(a, c);
        
        if(b > d)
            std::swap(b, d);
        
        if(b > c)
            std::swap(b, c);
    }
    else
    {
        if(a < b)
            std::swap(a, b);
        
        if(c < d)
            std::swap(c, d);
        
        if(a < c)
            std::swap(a, c);
        
        if(b < d)
            std::swap(b, d);
        
        if(b < c)
            std::swap(b, c);
    }
}

template void im::core_sort_4(int &a, int &b, int &c, int &d, SortDirection dir);
template void im::core_sort_4(float &a, float &b, float &c, float &d, SortDirection dir);
template void im::core_sort_4(double &a, double &b, double &c, double &d, SortDirection dir);

// Ensure that an angle lies in the range -PI < angle <= PI
template <typename T> T im::core_cyclic_modulus(T angle)
{
    if(angle>-CONST_PI && angle<=CONST_PI)
        return angle;
    
    if(angle<0.0)
        return angle - 2.0*CONST_PI * (double)( ((int)(angle / CONST_PI) - 1)/2 );
    else
        return angle - 2.0*CONST_PI * (double)( ((int)(angle / CONST_PI) + 1)/2 );
}

template float im::core_cyclic_modulus(float angle);
template double im::core_cyclic_modulus(double angle);

// Calculate log to the base of 2. Returns -1 if value is not a power of 2.
int im::core_integer_log2(uint64_t value)
{
    if(value==0)
        return -1;
    
    int logval = 0;
    uint64_t test = 1;
    
    while(test < value)
    {
        test *= 2;
        logval++;
    }
    
    if(value != (((uint64_t)1)<<logval))
        return -1;
    
    return logval;
}

template <typename TT> TT im::core_hypot(TT x, TT y)
{
    if(x<0)
        x = -x;
    if(y<0)
        y = -y;
    
    if(x>y)
    {
        TT r = y/x;
        return x * std::sqrt(1 + r*r);
    }
    else if(y>0)
    {
        TT r = x/y;
        return y * std::sqrt(1 + r*r);
    }
    else
        return x;
}

template float im::core_hypot(float x, float y);
template double im::core_hypot(double x, double y);
