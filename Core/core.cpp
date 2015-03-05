//
//  core.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

bool im::MtxRect::is_intersecting(MtxRect const &other) const
{
    int r0 = std::max(origin.row, other.origin.row);
    int r1 = std::min(origin.row + size.rows, other.origin.row + other.size.rows);
    if(r1 <= r0)
        return false;
    
    int c0 = std::max(origin.col, other.origin.col);
    int c1 = std::min(origin.col + size.cols, other.origin.col + other.size.cols);
    if(c1 <= c0)
        return false;
    
    return true;
}

bool im::MtxRect::is_within(MtxRect const &other) const
{
    if(origin.row < other.origin.row || origin.row + size.rows > other.origin.row + other.size.rows)
        return false;
    if(origin.col < other.origin.col || origin.col + size.cols > other.origin.col + other.size.cols)
        return false;
    
    return true;
}

im::MtxRect im::MtxRect::clip(MtxRect const &other) const
{
    int r0 = std::max(origin.row, other.origin.row);
    int r1 = std::min(origin.row + size.rows, other.origin.row + other.size.rows);
    int nrows = r1 - r0;
    int c0 = std::max(origin.col, other.origin.col);
    int c1 = std::min(origin.col + size.cols, other.origin.col + other.size.cols);
    int ncols = c1 - c0;
    
    if(nrows<1 || ncols<1)
        return MtxRect(origin, MtxSize(0,0));
    else
        return MtxRect(r0, c0, nrows, ncols);
}

void im::core_print_value(FILE *fp, uint8_t const &p)
{
    fprintf(fp, "%d", p);
}

void im::core_print_value(FILE *fp, int16_t const &p)
{
    fprintf(fp, "%d", p);
}

void im::core_print_value(FILE *fp, int const &p)
{
    fprintf(fp, "%d", p);
}

void im::core_print_value(FILE *fp, float const &p)
{
    fprintf(fp, "%g", p);
}

void im::core_print_value(FILE *fp, double const &p)
{
    fprintf(fp, "%g", p);
}

void im::core_print_value(FILE *fp, im::Cf const &p)
{
    if(p.imag()==0)
        fprintf(fp, "%g", p.real());
    else if(p.real()==0)
        fprintf(fp, "%gi", p.imag());
    else if(p.imag()>0)
        fprintf(fp, "%g+%gi", p.real(), p.imag());
    else
        fprintf(fp, "%g-%gi", p.real(), -p.imag());
}

void im::core_print_value(FILE *fp, im::Cd const &p)
{
    if(p.imag()==0)
        fprintf(fp, "%g", p.real());
    else if(p.real()==0)
        fprintf(fp, "%gi", p.imag());
    else if(p.imag()>0)
        fprintf(fp, "%g+%gi", p.real(), p.imag());
    else
        fprintf(fp, "%g-%gi", p.real(), -p.imag());
}

void im::core_write_value(FILE *fp, uint8_t const &p, char delimiter)
{
    fprintf(fp, "%d", p);
}

void im::core_write_value(FILE *fp, int16_t const &p, char delimiter)
{
    fprintf(fp, "%d", p);
}

void im::core_write_value(FILE *fp, int const &p, char delimiter)
{
    fprintf(fp, "%d", p);
}

void im::core_write_value(FILE *fp, float const &p, char delimiter)
{
    fprintf(fp, "%g", p);
}

void im::core_write_value(FILE *fp, double const &p, char delimiter)
{
    fprintf(fp, "%g", p);
}

void im::core_write_value(FILE *fp, im::Cf const &p, char delimiter)
{
    fprintf(fp, "%g%c%g", p.real(), delimiter, p.imag());
}

void im::core_write_value(FILE *fp, im::Cd const &p, char delimiter)
{
    fprintf(fp, "%g%c%g", p.real(), delimiter, p.imag());
}





