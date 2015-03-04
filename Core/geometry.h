//
//  geometry.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_geometry_h
#define ImpressiveAI_geometry_h

namespace im
{
    // Returns the coefficients of a line equation that passes through the thow points given
    // The coefficients vector has the form (a b c), for the line equation ax + by + c = 0
    template <typename TT> void core_line_through_points(VecView<TT> vcoefs, TT x1, TT y1, TT x2, TT y2);
    
    // Make similarity transform for homogeneous 2D transform
    template <typename TT> void core_make_3x3_similarity(MtxView<TT> m, TT scale, TT angle, TT tx, TT ty);
    
    // Make 3x3 skew-symmetric matrix
    template <typename TT> void core_make_3x3_skew_symmetric(MtxView<TT> m, TT dx, TT dy, TT dz);
    
    // Make 3x3 matrices for 3D transforms
    template <typename TT> void core_make_3x3_rotation_about_x(MtxView<TT> m, TT angle);
    template <typename TT> void core_make_3x3_rotation_about_y(MtxView<TT> m, TT angle);
    template <typename TT> void core_make_3x3_rotation_about_z(MtxView<TT> m, TT angle);
    template <typename TT> void core_make_3x3_rotation_axis_angle(MtxView<TT> m, VecView<TT> const &vaxis, TT angle);
    template <typename TT> void core_make_3x3_rotation_euler(MtxView<TT> m, TT roll, TT pitch, TT yaw);
    template <typename TT> void core_make_3x3_scale(MtxView<TT> m, TT sx, TT sy, TT sz);
    
    // Make 4x4 matrices for homogeneous 3D transforms
    template <typename TT> void core_make_4x4_rotation_about_x(MtxView<TT> m, TT angle);
    template <typename TT> void core_make_4x4_rotation_about_y(MtxView<TT> m, TT angle);
    template <typename TT> void core_make_4x4_rotation_about_z(MtxView<TT> m, TT angle);
    template <typename TT> void core_make_4x4_rotation_axis_angle(MtxView<TT> m, VecView<TT> const &vaxis, TT angle);
    template <typename TT> void core_make_4x4_rotation_euler(MtxView<TT> m, TT roll, TT pitch, TT yaw);
    template <typename TT> void core_make_4x4_scale(MtxView<TT> m, TT sx, TT sy, TT sz);
    template <typename TT> void core_make_4x4_translate(MtxView<TT> m, TT tx, TT ty, TT tz);
    
    
}

#endif
