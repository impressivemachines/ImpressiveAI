//
//  types.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/29/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_types_h
#define Metaphor_types_h

namespace im
{
    // Type identifiers for values used in the library
    enum CoreType {
        CoreTypeUndefined = 0,
        CoreTypeByte,
        CoreTypeShort,
        CoreTypeInt,
        CoreTypeFloat,
        CoreTypeDouble,
        CoreTypeComplexFloat,
        CoreTypeComplexDouble,
    };
    
    // Properties of common types
    
    // Usage:
    // TypeProperties<TT>::core_type_id - id number for type
    // TypeProperties<TT>::is_integer - true for integer types
    // TypeProperties<TT>::is_signed - true for signed types
    // TypeProperties<TT>::is_complex - true for complex number types
    // TypeProperties<TT>::Scalar - the inner type, e.g. float for std::complex<float>
    // TypeProperties<TT>::Promoted - becomes int for integer types with smaller range than int, else same as type
    // TypeProperties<TT>::Precise - the most accurate type of its class, e.g. short -> int, float -> double, std::complex<float> -> std::complex<double>
    // TypeProperties<TT>::epsilon() - smallest difference from 1
    // TypeProperties<TT>::highest() - largest possible positive number representable
    // TypeProperties<TT>::lowest() - most negative number representable
    
    template <typename TT>
    struct TypePropertiesBase
    {
        enum {
            core_type_id = CoreTypeUndefined,
            is_integer = std::numeric_limits<TT>::is_integer,
            is_signed = std::numeric_limits<TT>::is_signed,
            is_complex = false
        };
        
        typedef TT Scalar;
        typedef TT Promoted;
        typedef TT Precise;
        
        static inline Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon(); }
        static inline Scalar highest() { return (std::numeric_limits<Scalar>::max)(); }
        static inline Scalar lowest()  { return is_integer ? (std::numeric_limits<Scalar>::min)() : (-(std::numeric_limits<Scalar>::max)()); }
    };
    
    // Fallback
    template<typename TT> struct TypeProperties : TypePropertiesBase<TT>
    {};
    
    // CoreTypeByte
    template<> struct TypeProperties<uint8_t> : TypePropertiesBase<uint8_t>
    {
        enum { core_type_id = CoreTypeByte };
        typedef int Promoted;
        typedef int Precise;
    };
    
    // CoreTypeShort
    template<> struct TypeProperties<int16_t> : TypePropertiesBase<int16_t>
    {
        enum { core_type_id = CoreTypeShort };
        typedef int Promoted;
        typedef int Precise;
    };
    
    // CoreTypeInt
    template<> struct TypeProperties<int> : TypePropertiesBase<int>
    {
        enum { core_type_id = CoreTypeInt };
    };
    
    // CoreTypeFloat
    template<> struct TypeProperties<float> : TypePropertiesBase<float>
    {
        enum { core_type_id = CoreTypeFloat };
        
        typedef double Precise;
    };
    
    // CoreTypeDouble
    template<> struct TypeProperties<double> : TypePropertiesBase<double>
    {
        enum { core_type_id = CoreTypeDouble };
    };
    
    // CoreTypeComplexFloat
    template<> struct TypeProperties< std::complex<float> > : TypePropertiesBase< std::complex<float> >
    {
        enum {
            core_type_id = CoreTypeComplexFloat,
            is_integer = false,
            is_signed = true,
            is_complex = true
        };
        
        typedef float Scalar;
        typedef std::complex<double> Precise;
        
        static inline Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon(); }
        static inline Scalar highest() { return std::numeric_limits<Scalar>::max(); }
        static inline Scalar lowest()  { return -std::numeric_limits<Scalar>::max(); }
    };
    
    // CoreTypeComplexDouble
    template<> struct TypeProperties< std::complex<double> > : TypePropertiesBase< std::complex<double> >
    {
        enum {
            core_type_id = CoreTypeComplexDouble,
            is_integer = false,
            is_signed =  true,
            is_complex = true
        };
        
        typedef double Scalar;
        
        static inline Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon(); }
        static inline Scalar highest() { return std::numeric_limits<Scalar>::max(); }
        static inline Scalar lowest()  { return -std::numeric_limits<Scalar>::max(); }
    };
    
}

#endif
