//
//  quaternion.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_quaternion_h
#define Metaphor_quaternion_h

namespace im
{
    template <typename TT>
    class Quat
    {
    public:
        Quat() {}
        Quat(TT vx, TT vy, TT vz, TT vw) : x(vx), y(vy), z(vz), w(vw) {}
        Quat(VecView<TT> const &v3, TT vw) : x(v3(0)), y(v3(1)), z(v3(2)), w(vw) {}
        Quat(VecView<TT> const &v4) : x(v4(0)), y(v4(1)), z(v4(2)), w(v4(3)) {}
        
        void set_identity() { x = y = z = (TT)0; w = (TT)1; }
        
        // interpolate between two quaternions
        void from_iterpolation(Quat const &q1, Quat const &q2, TT alpha);
        
        // make a quaternion expressing the rotation between two 3-vectors
        void from_vector_rotation(VecView<TT> const &vorigin, VecView<TT> const &vdest);
        
        // create from euler angles
        void from_euler_angles(TT pitch, TT yaw, TT roll);
        
        // create from axis 3-vector and angle of rotation
        void from_axis_angle(VecView<TT> const &vaxis, TT angle);
        
        // create from 3x3 rotation matrix
        void from_rotation_matrix(MtxView<TT> const &mrot);
        
        // convert to euler angles
        void to_euler_angles(TT &pitch, TT &yaw, TT &roll);
        
        // convert to 3-vector axis of rotation and angle
        void to_axis_angle(VecView<TT> vaxis, TT &angle);
        
        // convert to 3x3 rotation matrix
        void to_rotation_matrix(MtxView<TT> mrot);
        
        TT magnitude() const { return std::sqrt(magnitude_squared()); }
        TT magnitude_squared() const { return x*x + y*y + z*z + w*w; }
        
        // normalize
        Quat unit() const { return (*this)/magnitude(); }
        
        // conjugate
        Quat conj() const { return Quat(-x, -y, -z, w); }
        
        // inverse of the quaternion
        Quat inverse() const { return conj() / magnitude_squared(); }
        
        TT dot_product(Quat const &q) { return x * q.x + y * q.y + z * q.z + w * q.w; }
        
        // operators
        Quat operator+() const { return *this; }
        Quat operator-() const { return Quat(-x, -y, -z, -w); }
        Quat operator+(Quat const &q) { return Quat(x+q.x, y+q.y, z+q.z, w+q.w); }
        Quat operator-(Quat const &q) { return Quat(x-q.x, y-q.y, z-q.z, w-q.w); }
        Quat &operator+=(Quat const &q) { x+=q.x; y+=q.y; z+=q.z; w+=q.w; return *this; }
        Quat &operator-=(Quat const &q) { x-=q.x; y-=q.y; z-=q.z; w-=q.w; return *this; }
        
        // combine solid body rotations
        Quat operator*(Quat const &q);
        Quat operator*=(Quat const &q) { *this = *this * q; return *this; }
        
        // scaling
        Quat operator*(TT f) const { return Quat(x*f, y*f, z*f, w*f); }
        Quat operator/(TT f) const { return Quat(x/f, y/f, z/f, w/f); }
        Quat &operator*=(TT f) { x*=f; y*=f; z*=f; w*=f; return *this; }
        Quat &operator/=(TT f) { x/=f; y/=f; z/=f; w/=f; return *this; }
        
        // compare
        bool operator==(Quat const &q) const { return x==q.x && y==q.y && z==q.z && w==q.w; }
        bool operator!=(Quat const &q) const { return !operator==(q); }
        
    public:
        TT x, y, z, w;
    };

    template <typename TT> Quat<TT> operator*(const TT &f, const Quat<TT> &q)
    {
        return Quat<TT>(q.x*f, q.y*f, q.z*f, q.w*f);
    }
}

#endif
