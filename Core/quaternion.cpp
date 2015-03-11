//
//  quaternion.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

// interpolate between two quaternions
template <typename TT>
void im::Quat<TT>::from_iterpolation(Quat const &qstart, Quat const &qend, TT alpha)
{
    alpha = std::min((TT)1, alpha);
    alpha = std::max((TT)0, alpha);
    
    Quat<TT> q1 = qstart.unit();
    Quat<TT> q2 = qend.unit();
    
    if((q1 - q2).magnitude_squared() > (q1 + q2).magnitude_squared())
        q2 = -q2;
    
    // omega is half the rotation angle
    TT cos_omega = std::min((TT)1, q1.dot_product(q2));
    TT sin_omega = std::sqrt((TT)1 - cos_omega * cos_omega); // always positive
    
    if(sin_omega<(TT)0.0001)
    {
        if(cos_omega>(TT)0)
        {
            // angle close to zero so just lerp
            *this = q1 * ((TT)1 - alpha) + q2 * alpha;
        }
        else
        {
            // angle close to pi so ill defined
            x = y = z = w = (TT)0;
        }
    }
    else
    {
        // slerp
        TT omega = std::acos(cos_omega);
        *this = q1 * (std::sin(((TT)1 - alpha) * omega)/sin_omega) + q2 * (std::sin(alpha * omega)/sin_omega);
    }
}

// make a quaternion expressing the rotation between two 3-vectors
template <typename TT>
void im::Quat<TT>::from_vector_rotation(VecView<TT> const &vorigin, VecView<TT> const &vdest)
{
    IM_CHECK_VALID(vorigin);
    IM_CHECK_VECTOR_SIZE(vorigin, 3);
    IM_CHECK_VALID(vdest);
    IM_CHECK_VECTOR_SIZE(vdest, 3);
    
    TT vc[3];
    vc[0] = vorigin(1) * vdest(2) - vorigin(2) * vdest(1);
    vc[1] = vorigin(2) * vdest(0) - vorigin(0) * vdest(2);
    vc[2] = vorigin(0) * vdest(1) - vorigin(1) * vdest(0);

    TT mag = std::sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);

    if(mag < TypeProperties<TT>::epsilon())
    {
        x = (TT)0;
        y = (TT)0;
        z = (TT)0;
        w = (TT)1;
    }
    else
    {
        TT c2 = vorigin(0) * vdest(0) + vorigin(1) * vdest(1) + vorigin(2) * vdest(2);
        c2 /= std::sqrt(vorigin(0) * vorigin(0) + vorigin(1) * vorigin(1) + vorigin(2) * vorigin(2));
        c2 /= std::sqrt(vdest(0) * vdest(0) + vdest(1) * vdest(1) + vdest(2) * vdest(2));
        
        TT fs  = std::sqrt(((TT)1 - c2)/(TT)2);
        TT fc  = std::sqrt(((TT)1 + c2)/(TT)2);
        
        x = vc[0] * fs / mag;
        y = vc[1] * fs / mag;
        z = vc[2] * fs / mag;
        w = fc;
    }
}

// create from euler angles
template <typename TT>
void im::Quat<TT>::from_euler_angles(TT pitch, TT yaw, TT roll)
{
    // This is using the Tait-Bryan angles with composition order Y_yaw X_pitch Z_roll
    
    Quat<TT> qx, qy, qz;
    TT u[3];
    
    u[0] = 1;
    u[1] = 0;
    u[2] = 0;
    
    qx.from_axis_angle(VecView<TT>(3, 1, u), pitch); // about x
    
    u[0] = 0;
    u[1] = 1;
    u[2] = 0;
    
    qy.from_axis_angle(VecView<TT>(3, 1, u), yaw); // about y
    
    u[0] = 0;
    u[1] = 0;
    u[2] = 1;
    
    qz.from_axis_angle(VecView<TT>(3, 1, u), roll); // about z
    
    *this = qy * (qx * qz);
}

// create from axis 3-vector and angle of rotation
template <typename TT>
void im::Quat<TT>::from_axis_angle(VecView<TT> const &vaxis, TT angle)
{
    IM_CHECK_VALID(vaxis);
    IM_CHECK_VECTOR_SIZE(vaxis, 3);
    
    angle *= (TT)0.5;
    
    TT fs = std::sin(angle);
    TT fc = std::cos(angle);
    
    TT magsq = vaxis(0) * vaxis(0) + vaxis(1) * vaxis(1) + vaxis(2) * vaxis(2);
    if(magsq==(TT)0)
    {
        x = y = z = (TT)0;
    }
    else
    {
        TT mag = std::sqrt(magsq);
        x = fs * vaxis(0)/mag;
        y = fs * vaxis(1)/mag;
        z = fs * vaxis(2)/mag;
    }
    w = fc;
}

// create from 3x3 rotation matrix
template <typename TT>
void im::Quat<TT>::from_rotation_matrix(MtxView<TT> const &mrot)
{
    // See "Quaternions" by Ken Shoemake, U. Penn.
    
    IM_CHECK_VALID(mrot);
    IM_CHECK_MATRIX_SIZE(mrot, 3, 3);

    TT trace = mrot(0,0) + mrot(1,1) + mrot(2,2);
    
    if(trace>(TT)0)
    {
        TT root = std::sqrt(trace +  (TT)1);
        w = root * (TT)0.5;
        root = (TT)0.5/root;
        
        x = (mrot(2,1) - mrot(1,2)) * root;
        y = (mrot(0,2) - mrot(2,0)) * root;
        z = (mrot(1,0) - mrot(0,1)) * root;
    }
    else
    {
        static int const cycle[3] = { 1, 2, 0 };
        
        int i;
        i = mrot(1,1) > mrot(0,0) ? 1 : 0;
        i = mrot(2,2) > mrot(i,i) ? 2 : i;
        
        int j = cycle[i];
        int k = cycle[j];
        
        TT root = std::sqrt(mrot(i,i) - mrot(j,j) - mrot(k,k) + (TT)1);
        
        TT a = (TT)0.5 * root;
        if(i==0)
            x = a;
        else if(i==1)
            y = a;
        else
            z = a;
        
        if(root != (TT)0)
            root = (TT)0.5 / root;

        a = (mrot(i,j) + mrot(j,i)) * root;
        if(j==0)
            x = a;
        else if(j==1)
            y = a;
        else
            z = a;
        
        a = (mrot(i,k) + mrot(k,i)) * root;
        if(k==0)
            x = a;
        else if(k==1)
            y = a;
        else
            z = a;
        
        w = (mrot(k,j) - mrot(j,k)) * root;
    }
}

// convert to euler angles
template <typename TT>
void im::Quat<TT>::to_euler_angles(TT &pitch, TT &yaw, TT &roll)
{
    TT u[9];
    MtxView<TT> mR(3, 3, 3, 1, u);
    to_rotation_matrix(mR);
    
    // This is using the Tait-Bryan angles with composition order Y_yaw X_pitch Z_roll
    // given r, p, y, matrix is:
    // [ cy*cr+sy*sr*sp      -cy*sr+sy*cr*sp     sy*cp;
    //   sr*cp               cr*cp               -sp;
    //   -sy*cr+cy*sr*sp     sy*sr+cy*cr*sp      cy*cp ]
    
    pitch = std::asin(-mR(1,2));
    if((mR(0,2)!=(TT)0 || mR(2,2)!=(TT)0)
       && (mR(1,0)!=(TT)0 || mR(1,1)!=(TT)0))
    {
        yaw = std::atan2(mR(0,2), mR(2,2));
        roll = std::atan2(mR(1,0), mR(1,1));
    }
    else
    {
        roll = (TT)0;
        yaw = std::atan2(-mR(2,0), mR(0,0));
    }
}

// convert to 3-vector axis of rotation and angle
template <typename TT>
void im::Quat<TT>::to_axis_angle(VecView<TT> vaxis, TT &angle)
{
    IM_CHECK_VALID(vaxis);
    IM_CHECK_VECTOR_SIZE(vaxis, 3);
    
    Quat<TT> q = unit();
    angle = std::acos(q.w);
    
    if(std::abs(angle)<=TypeProperties<TT>::epsilon())
    {
        // default direction to keep axis normalized
        vaxis(0) = (TT)1;
        vaxis(1) = (TT)0;
        vaxis(2) = (TT)0;
    }
    else
    {
        TT fs = std::sin(angle);
        
        vaxis(0) = q.x / fs;
        vaxis(1) = q.y / fs;
        vaxis(2) = q.z / fs;
        
        TT mag = std::sqrt(vaxis(0) * vaxis(0) + vaxis(1) * vaxis(1) + vaxis(2) * vaxis(2));
        
        vaxis(0) /= mag;
        vaxis(1) /= mag;
        vaxis(2) /= mag;
    }
    
    angle *= (TT)2;
}

// convert to 3x3 rotation matrix
template <typename TT>
void im::Quat<TT>::to_rotation_matrix(MtxView<TT> mrot)
{
    // See "Quaternions" by Ken Shoemake, U. Penn.
    
    IM_CHECK_VALID(mrot);
    IM_CHECK_MATRIX_SIZE(mrot, 3, 3);
    
    TT wx = (TT)2 * w * x;
    TT wy = (TT)2 * w * y;
    TT wz = (TT)2 * w * z;
    TT xx = (TT)2 * x * x;
    TT xy = (TT)2 * x * y;
    TT xz = (TT)2 * x * z;
    TT yy = (TT)2 * y * y;
    TT yz = (TT)2 * y * z;
    TT zz = (TT)2 * z * z;
    
    mrot(0,0) = (TT)1 - yy - zz;
    mrot(0,1) = xy - wz;
    mrot(0,2) = xz + wy;
    
    mrot(1,0) = xy + wz;
    mrot(1,1) = (TT)1 - xx - zz;
    mrot(1,2) = yz - wx;
    
    mrot(2,0) = xz - wy;
    mrot(2,1) = yz + wx;
    mrot(2,2) = (TT)1 - xx - yy;
}

template <typename TT>
im::Quat<TT> im::Quat<TT>::operator*(Quat const &q)
{
    // v = (v cross qv) + qv * w + v * qw
    // w = w * qw - (v dot qv)
    
    TT newx = y * q.z - z * q.y + x * q.w + w * q.x;
    TT newy = z * q.x - x * q.z + y * q.w + w * q.y;
    TT newz = x * q.y - y * q.x + z * q.w + w * q.z;
    TT neww = w * q.w - x * q.x - y * q.y - z * q.z;
    
    return Quat<TT>(newx, newy, newz, neww);
}

template class im::Quat<float>;
template class im::Quat<double>;

