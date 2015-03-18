//
//  meta_core.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/24/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_meta_core_h
#define Metaphor_meta_core_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <stdexcept>
#include <chrono>
#include <memory>
#include <cmath>
#include <vector>
#include <random>
#include <complex>

// Math constants to double precision

#define CONST_E             2.71828182845904523536
#define CONST_2PI           6.28318530717958647692
#define CONST_PI            3.14159265358979323846
#define CONST_PI_2          1.57079632679489661923
#define CONST_PI_4          0.78539816339744830961
#define CONST_PI_180        0.01745329251994329576
#define CONST_G             1.61803398874989484820
#define CONST_LOG10E        0.43429448190325182765
#define CONST_LOG2E         1.44269504088896340735
#define CONST_LN10          2.30258509299404568401
#define CONST_LN2           0.69314718055994530941
#define CONST_1_PI          0.31830988618379067153
#define CONST_2_PI          0.63661977236758134307
#define CONST_180_PI        57.2957795130823208767
#define CONST_2_SQRTPI      1.12837916709551257389
#define CONST_SQRT2         1.41421356237309504880
#define CONST_SQRT1_2       0.70710678118654752440
#define CONST_SQRT3         1.73205080756887729352
#define CONST_SQRT3_2       0.86602540378443864676
#define CONST_SQRT2PI       2.50662827463100050241

// Error checking
#define IM_THROW_GENERIC                    throw std::logic_error("Generic assert failed.")
#define IM_THROW_NULL                       throw std::logic_error("Attempt to access NULL pointer.")
#define IM_THROW_NO_INIT                    throw std::logic_error("Uninitialized object or data.")
#define IM_THROW_BOUNDS                     throw std::out_of_range("Bounds check failed.")
#define IM_THROW_MATRIX                     throw std::domain_error("Incorrect matrix format.")
#define IM_THROW_ARGUMENTS                  throw std::invalid_argument("Invalid arguments.")
#define IM_THROW_FILE_ERROR                 throw std::runtime_error("File error.")
#define IM_THROW_FILE_FORMAT_ERROR          throw std::runtime_error("File format error.")
#define IM_THROW_NO_SOLUTION                throw std::runtime_error("No solution found.")
#define IM_THROW_BAD_DATA                   throw std::runtime_error("Unexpected data values encountered.")
#define IM_THROW_END_OF_FILE                throw std::runtime_error("End of file encountered.")
#define IM_THROW_ALLOC                      throw std::runtime_error("Memory allocation failed.")
#define IM_THROW_NO_IMPL                    throw std::runtime_error("Not implemented.")

#define IM_CHECK(e)                         do { if(!(e)) IM_THROW_GENERIC; } while(false)
#define IM_CHECK_ARGS(e)                    do { if(!(e)) IM_THROW_ARGUMENTS; } while(false)
#define IM_CHECK_NULL(p)                    do { if(p==NULL) IM_THROW_NULL; } while(false)

// check that i >= lower && i < upper
#define IM_CHECK_BOUNDS(i, lower, upper)    do { if(!((i)>=(lower) && (i)<(upper))) IM_THROW_BOUNDS; } while(false)
// check that i >= lower
#define IM_CHECK_LOWER_BOUNDS(i, lower)      do { if(!((i)>=(lower))) IM_THROW_BOUNDS; } while(false)
// check that i < upper
#define IM_CHECK_UPPER_BOUNDS(i, upper)      do { if(!((i)<(upper))) IM_THROW_BOUNDS; } while(false)
// check object valid
#define IM_CHECK_VALID(m)                   do { if(!((m).is_valid())) IM_THROW_MATRIX; } while(false)
// check matrix size
#define IM_CHECK_MATRIX_SIZE(m,r,c)          do { if(!((m).rows()==(r) && (m).cols()==(c))) IM_THROW_MATRIX; } while(false)
// check vector size
#define IM_CHECK_VECTOR_SIZE(m,r)           do { if(!((m).rows()==(r))) IM_THROW_MATRIX; } while(false)
// check matrix is square
#define IM_CHECK_MATRIX_SQUARE(m)            do { if((m).cols()!=(m).rows()) IM_THROW_MATRIX; } while(false)
// check two matrices are same size
#define IM_CHECK_MATRIX_SIZES_MATCH(m1,m2)   do { if(!((m1).rows()==(m2).rows() && (m1).cols()==(m2).cols())) IM_THROW_MATRIX; } while(false)
// check two vectors are same size
#define IM_CHECK_VECTOR_SIZES_MATCH(m1,m2)   do { if(!((m1).rows()==(m2).rows())) IM_THROW_MATRIX; } while(false)

// Time critical debug-only checks for inner loops or element get/set
#ifdef DEBUG

// check that i >= lower && i < upper
#define IM_DEBUG_ONLY_CHECK_BOUNDS(i, lower, upper)     do { if(!((i)>=(lower) && (i)<(upper))) IM_THROW_BOUNDS; } while(false)
// check that i >= lower
#define IM_DEBUG_ONLY_CHECK_LOWER_BOUNDS(i, lower)      do { if(!((i)>=(lower))) IM_THROW_BOUNDS; } while(false)
// check that i < upper
#define IM_DEBUG_ONLY_CHECK_UPPER_BOUNDS(i, upper)      do { if(!((i)<(upper))) IM_THROW_BOUNDS; } while(false)
// check NULL
#define IM_DEBUG_ONLY_CHECK_NULL(p)                     do { if(p==NULL) IM_THROW_GENERIC; } while(false)

#else

#define IM_DEBUG_ONLY_CHECK_BOUNDS(i, lower, upper)
#define IM_DEBUG_ONLY_CHECK_LOWER_BOUNDS(i, lower)
#define IM_DEBUG_ONLY_CHECK_UPPER_BOUNDS(i, upper)
#define IM_DEBUG_ONLY_CHECK_NULL(p)

#endif

// typical alignment for 128 bit (four float) SSE words

// k must be a power of 2
#define IM_SSE_ALIGN                    16
#define IM_IS_ALIGNED_K(p,k)            (( ((intptr_t)(p)) & (k-1))==0)
#define IM_IS_ALIGNED_SSE(p)            IM_IS_ALIGNED_K(p,IM_SSE_ALIGN)
#define IM_IS_SIZED_SSE(c)              (((c) & (IM_SSE_ALIGN-1))==0)

namespace im
{
    // Complex types
    typedef std::complex<float> Cf;
    typedef std::complex<double> Cd;
    
    // Complex number helpers
    
    // Real part of the value as scalar type
    template <typename T> inline T core_real(T const &s) { return s; }
    template <typename T> inline T core_real(std::complex<T> const &c) { return c.real(); }
    
    // Imaginary part of the value as scalar type
    template <typename T> inline T core_imag(T const &s) { return (T)0.0; }
    template <typename T> inline T core_imag(std::complex<T> const &c) { return c.imag(); }
    
    // Sum of absolute values |re|+|im| (not magnitude)
    template <typename T> inline T core_sumabs(T const &s) { return std::abs(s); }
    template <typename T> inline T core_sumabs(std::complex<T> const &c) { return std::abs(c.real()) + std::abs(c.imag()); }
    
    // Conjugate if complex
    template <typename T> inline T core_conj(T const &s) { return s; }
    template <typename T> inline std::complex<T> core_conj(std::complex<T> const &c) { return std::conj(c); }
    
    // Return 1 if val>0, 0 if val==0, and -1 if val<0
    template <typename T> int core_sgn(T const &val) { return (val > T(0)) - (val < T(0)); }
    
    // Forces x to be truly represented by the type, e.g. double gets to be 64 bits, without extended precision
    template <typename T> T core_coerce(const T x)
    {
        volatile T y;
        y = x;
        return y;
    }
    
    template <typename T> int core_round(const T x) { return (int)std::floor(x+(T)0.5); }
    
    // Structure for specifying matrix element locations
    struct MtxLoc
    {
        MtxLoc() {}
        MtxLoc(int r, int c) : row(r), col(c) {}
        
        int row, col;
    };
    
    struct MtxSize
    {
        MtxSize() {}
        MtxSize(int nrows, int ncols) : rows(nrows), cols(ncols) {}
        
        bool is_empty() const { return rows<1 || cols<1; }
        MtxSize clip(MtxSize const &other) const { return MtxSize(std::min(rows, other.rows), std::min(cols, other.cols)); }
        void normalize() { rows = std::abs(rows); cols = std::abs(cols); }
        
        int rows, cols;
    };
    
    struct MtxRect
    {
        MtxRect() {}
        MtxRect(int r, int c, int nrows, int ncols) : origin(r, c), size(nrows, ncols) {}
        MtxRect(MtxLoc org, MtxSize sz) : origin(org), size(sz) {}
        
        bool is_empty() const { return size.is_empty(); }
        bool is_intersecting(MtxRect const &other) const;
        bool is_within(MtxRect const &other) const;
        MtxRect clip(MtxRect const &other) const;
        
        MtxLoc origin;
        MtxSize size;
    };
    
    
    
    
    // For printing, typically in matlab format
    template <typename TT> void core_print_value(FILE *fp, TT const &p)
    {
        fprintf(fp, "??"); // unknown type fallback case
    }
    
    // Supported types - see core.cpp
    void core_print_value(FILE *fp, uint8_t const &p);
    void core_print_value(FILE *fp, int16_t const &p);
    void core_print_value(FILE *fp, int const &p);
    void core_print_value(FILE *fp, float const &p);
    void core_print_value(FILE *fp, double const &p);
    void core_print_value(FILE *fp, im::Cf const &p);
    void core_print_value(FILE *fp, im::Cd const &p);
    
    // For printing to CSV files
    template <typename TT> void core_write_value(FILE *fp, TT const &p, char delimiter)
    {
        fprintf(fp, "??"); // unknown type fallback case
    }
    
    // Supported types - see core.cpp
    void core_write_value(FILE *fp, uint8_t const &p, char delimiter);
    void core_write_value(FILE *fp, int16_t const  &p, char delimiter);
    void core_write_value(FILE *fp, int const &p, char delimiter);
    void core_write_value(FILE *fp, float const &p, char delimiter);
    void core_write_value(FILE *fp, double const &p, char delimiter);
    void core_write_value(FILE *fp, im::Cf const &p, char delimiter);
    void core_write_value(FILE *fp, im::Cd const &p, char delimiter);
    
    enum SortDirection
    {
        SortDirectionAscending, SortDirectionDescending
    };
    
    enum PadMode
    {
        PadModeZero, // pad with zero
        PadModeExtend, // pad with repeat of edge element
        PadModeWrap, // pad by wrapping from other side
        PadModeReflect // pad by reflecting back into interior
    };
    
    enum TriMode
    {
        TriMode_L,          // lower triangular
        TriMode_U           // upper triangular
    };
    
    enum TransMode
    {
        TransMode_N,    // no traspose
        TransMode_T,    // transpose
        TransMode_C     // conjugate transpose
    };
    
    enum DiagMode
    {
        DiagMode_N,         // normal diagonal
        DiagMode_U          // unit diagonal
    };
    
    enum SideMode
    {
        SideMode_L,
        SideMode_R
    };
}

#include "types.h"
#include "mem.h"
#include "block_loop.h"
#include "block_generic.h"
#include "vector_view.h"
#include "matrix_view.h"
#include "block_blas1.h"
#include "block_blas2.h"
#include "block_blas3.h"
#include "block_unary.h"
#include "block_binary.h"
#include "block_reduce.h"
#include "block_convert.h"
#include "block_sort.h"
#include "block_tools.h"
#include "random.h"
#include "vector.h"
#include "matrix.h"
#include "convolve.h"
#include "utils.h"
#include "matrix_io.h"
#include "permute.h"
#include "decomp_tools.h"
#include "decomp_lu.h"
#include "decomp_qr.h"
#include "decomp_cholesky.h"
#include "decomp_tridiag.h"
#include "decomp_eigen.h"
#include "decomp_svd.h"
#include "solve.h"
#include "stats.h"
#include "fft.h"
#include "quaternion.h"
#include "geometry.h"
#include "opt_linemin.h"
#include "opt_levmarq.h"
#include "opt_powell.h"
#include "opt_conjgrad.h"
#include "opt_bfgs.h"
#include "opt_stochastic.h"

#endif
