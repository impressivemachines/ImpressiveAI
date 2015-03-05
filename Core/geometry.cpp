//
//  geometry.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::core_line_through_points(VecView<TT> vcoefs, TT x1, TT y1, TT x2, TT y2)
{
    IM_CHECK_VALID(vcoefs);
    IM_CHECK_VECTOR_SIZE(vcoefs, 3);
    
    TT a = y1-y2;
    TT b = x2-x1;
    TT c = x1*y2 - x2*y1;
    
    TT norm = std::sqrt(a*a+b*b);
    if(norm==(TT)0)
    {
        vcoefs(0) = 0;
        vcoefs(1) = -1;
        vcoefs(2) = y1;
    }
    else
    {
        vcoefs(0) = a / norm;
        vcoefs(1) = b / norm;
        vcoefs(2) = c / norm;
    }
}

#define INST(TT) template void im::core_line_through_points(VecView<TT> vcoefs, TT x1, TT y1, TT x2, TT y2)
INST(float); INST(double);
#undef INST

// Make similarity transform for 2D using a 3x3 homogeneous matrix
template <typename TT>
void im::core_make_3x3_similarity(MtxView<TT> m, TT scale, TT angle, TT tx, TT ty)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
}

// Make 3x3 skew-symmetric matrix
template <typename TT>
void im::core_make_3x3_skew_symmetric(MtxView<TT> m, TT dx, TT dy, TT dz)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    
    m(0,0) = (TT)0; m(0,1) = -dz;   m(0,2) = dy;
    m(1,0) = dz;    m(1,1) = (TT)0; m(1,2) = -dx;
    m(2,0) = -dy;   m(2,1) = dx;    m(2,2) = (TT)0;
}

// Make 3x3 rotation matrices for 3D transforms
template <typename TT>
void im::core_make_3x3_rotation_about_x(MtxView<TT> m, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    
    TT s = sin(angle);
    TT c = cos(angle);
    
    m = (TT)0;
    m.diag() = (TT)1;
    
    m(1,1) = c;
    m(1,2) = -s;
    m(2,1) = s;
    m(2,2) = c;
}

template <typename TT>
void im::core_make_3x3_rotation_about_y(MtxView<TT> m, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    
    TT s = sin(angle);
    TT c = cos(angle);
    
    m = (TT)0;
    m.diag() = (TT)1;
    
    m(0,0) = c;
    m(0,2) = s;
    m(2,0) = -s;
    m(2,2) = c;
}

template <typename TT>
void im::core_make_3x3_rotation_about_z(MtxView<TT> m, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    
    TT s = sin(angle);
    TT c = cos(angle);
    
    m = (TT)0;
    m.diag() = (TT)1;
    
    m(0,0) = c;
    m(0,1) = -s;
    m(1,0) = s;
    m(1,1) = c;
}

template <typename TT>
void im::core_make_3x3_rotation_axis_angle(MtxView<TT> m, VecView<TT> const &vaxis, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_VALID(vaxis);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    IM_CHECK_VECTOR_SIZE(vaxis, 3);
    
    TT s = sin(angle);
    TT c = cos(angle);
    TT mc = 1 - c;
    
    TT v[3];
    v[0] = vaxis(0);
    v[1] = vaxis(1);
    v[2] = vaxis(2);
    
    TT mag2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    
    if(mag2==0.0)
    {
        v[0] = (TT)1;
        v[1] = (TT)0;
        v[2] = (TT)0;
    }
    else
    {
        mag2 = std::sqrt(mag2);
        v[0] /= mag2;
        v[1] /= mag2;
        v[2] /= mag2;
    }
    
    m(0,0) = c + v[0] * v[0] * mc;
    m(0,1) = v[0] * v[1] * mc - v[2] * s;
    m(0,2) = v[0] * v[2] * mc + v[1] * s;
    
    m(1,0) = v[0] * v[1] * mc + v[2] * s;
    m(1,1) = c + v[1] * v[1] * mc;
    m(1,2) = v[1] * v[2] * mc - v[0] * s;
    
    m(2,0) = v[0] * v[2] * mc - v[1] * s;
    m(2,1) = v[1] * v[2] * mc + v[0] * s;
    m(2,2) = c + v[2] * v[2] * mc;
}

template <typename TT>
void im::core_make_3x3_rotation_euler(MtxView<TT> m, TT roll, TT pitch, TT yaw)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    
    TT ca = cos(yaw);
    TT sa = sin(yaw);
    TT cb = cos(pitch);
    TT sb = sin(pitch);
    TT cc = cos(roll);
    TT sc = sin(roll);
    
    m(0,0) = ca * cb;
    m(0,1) = ca * sb * sc - sa * cc;
    m(0,2) = ca * sb * cc + sa * sc;
    
    m(1,0) = sa * cb;
    m(1,1) = sa * sb * sc + ca * cc;
    m(1,2) = sa * sb * cc - ca * sc;
    
    m(2,0) = -sb;
    m(2,1) = cb * sc;
    m(2,2) = cb * cc;
}

template <typename TT>
void im::core_make_3x3_scale(MtxView<TT> m, TT sx, TT sy, TT sz)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 3, 3);
    
    m(0,0) = sx;
    m(0,1) = (TT)0;
    m(0,2) = (TT)0;
    
    m(1,0) = (TT)0;
    m(1,1) = sy;
    m(1,2) = (TT)0;
    
    m(2,0) = (TT)0;
    m(2,1) = (TT)0;
    m(2,2) = sz;
}

// Make 4x4 matrices for homogeneous 3D transforms
template <typename TT>
void im::core_make_4x4_rotation_about_x(MtxView<TT> m, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);
    core_make_3x3_rotation_about_x(m.block(0,0,3,3), angle);
    m.block(3,0,1,3) = (TT)0;
    m.block(0,3,3,1) = (TT)0;
    m(3,3) = (TT)1;
}

template <typename TT>
void im::core_make_4x4_rotation_about_y(MtxView<TT> m, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);
    core_make_3x3_rotation_about_y(m.block(0,0,3,3), angle);
    m.block(3,0,1,3) = (TT)0;
    m.block(0,3,3,1) = (TT)0;
    m(3,3) = (TT)1;
}

template <typename TT>
void im::core_make_4x4_rotation_about_z(MtxView<TT> m, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);
    core_make_3x3_rotation_about_z(m.block(0,0,3,3), angle);
    m.block(3,0,1,3) = (TT)0;
    m.block(0,3,3,1) = (TT)0;
    m(3,3) = (TT)1;
}

template <typename TT>
void im::core_make_4x4_rotation_axis_angle(MtxView<TT> m, VecView<TT> const &vaxis, TT angle)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);

    core_make_3x3_rotation_axis_angle(m.block(0,0,3,3), vaxis, angle);
    m.block(3,0,1,3) = (TT)0;
    m.block(0,3,3,1) = (TT)0;
    m(3,3) = (TT)1;
}

template <typename TT>
void im::core_make_4x4_rotation_euler(MtxView<TT> m, TT roll, TT pitch, TT yaw)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);
    
    core_make_3x3_rotation_euler(m.block(0,0,3,3), roll, pitch, yaw);
    m.block(3,0,1,3) = (TT)0;
    m.block(0,3,3,1) = (TT)0;
    m(3,3) = (TT)1;
}

template <typename TT>
void im::core_make_4x4_scale(MtxView<TT> m, TT sx, TT sy, TT sz)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);
    
    core_make_3x3_scale(m.block(0,0,3,3), sx, sy, sz);
    m.block(3,0,1,3) = (TT)0;
    m.block(0,3,3,1) = (TT)0;
    m(3,3) = (TT)1;
}

template <typename TT>
void im::core_make_4x4_translate(MtxView<TT> m, TT tx, TT ty, TT tz)
{
    IM_CHECK_VALID(m);
    IM_CHECK_MATRIX_SIZE(m, 4, 4);
    
    m.block(0,0,4,3) = (TT)0;
    m.diag() = (TT)1;
    m(0,3) = tx;
    m(1,3) = ty;
    m(2,3) = tz;
}
