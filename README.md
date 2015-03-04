# ImpressiveAI
Impressive Machines's C++ Library for Machine Learning and Computer Vision

The eventual aim of this library is to create a self-contained set of tools for many applications in machine learning and computer vision. It is broken into a number of sections which are useful in isolation. The main part of the library is "core" which is a library for linear algebra and numerical computation. This part currently has a lot of functionality. Other parts are under development and will be added in the future.

## Core
Core provides many common tools for linear algebra and statistics and is built on very flexible matrix and vector objects. The features include:

* Light-weight matrix and vector views which provide a unified, powerful, and simple way to refer to one and two dimensional memory in a multitude of configurations.
* A full suite of conforming BLAS operations supporting BLAS 1, 2 and 3 linear algebra manipulations. These make use of matrix and vector views for their arguments.
* Optimized low-level matrix block functions which work on vector and matrix types to implement common operations, such as math, data sorting, conversion, statistics and reduce operations.
* Matrix and vector memory objects which manage shared memory and can be used to express a multitude of elaborate mathematical operations neatly and efficiently.
* Real and complex types are supported for many math operations.
* Typical matrix decompositions:
    * LU decomposition.
    * Cholesky LDLT and LLT decomposition.
    * QR decomposition with and without pivoting.
    * Tridiagonal decomposition.
    * Eigenvalue decomposition for symmetric matrices.
    * Jacobi SVD with QR preconditioning.
    * Special case 2x2 Eigen and SVD decompositions.
* Fast convolution for 1D and 2D signals.
* Matrices and vectors can be saved and printed in a variety of formats, including Matlab text format.
* Convenient interfaces for random number generation.
* Solution to quadratic and cubic polynomials.
* Principal components analysis.
* Estimation of covariances.
* Generation of multidimensional Gaussian noise with specified covariance.
* Computation of multivariate Gaussians and Mahalanobis distance.
* Arbitrary radix FFT in 1D and 2D.
* Least squares and robust 2D line fitting.
* Multi-dimensional regression.
* Geometric 2D and 3D rotations and transforms.

Some features that are under development:
* Line minimization / root finding.
* Optimization of functions.
* Tensors.
* Quaternions.

The matrix and vector view objects have independent signed row and columns strides and so they can operate on many data formats and with any data type - they do not require classical row or column major memory organization. All operations support positive or negative row or column strides. Therefore they can refer to memory which is reversed, decimated, or fragmented, such as all red pixels of an RGB image with arbitrary row stride which is stored in bottom up order. 

Decimating, reversing, or extracting rows, columns, blocks, or diagonals of matrices is therefore a constant time operation. Most operations in the library involve passing only the small view structures as arguments to functions (often by reference) and unnecessary data copying is avoided. Matrix and vector classes exist to dynamically create memory blocks where needed and these make use of standard shared pointers. Alternatively, memory can be statically allocated by the user and trivially wrapped into matrix view objects for use with any library functions.

### Dependencies
There are no external dependencies for the core library.

### Philosophy
There are many linear algebra libraries available that are built using C++. Many of the popular ones make extensive use of template meta-programming to increase efficiency. This library takes a different approach. Not many developers really understand C++ templates well and so complex template code is hard to debug and maintain. 

This library aims towards simplicity of code design. It does this without giving up too much run time efficiency by making extensive use of light weight references to matrix memory. It avoids unneccesary data copies by using shared pointers where possible. The main goal was to make the code easy to understand and still quite fast. Although our library is a completely new code base, the primary reference source for creating many of the algorithms was the GNU scientific library and so the chosen routines are generally expected to be of good quality and have so far performed well in testing.

### Examples
    #include "imp_core.h"
    using namespace im;
    
    // Fill matrix matlab style
    Mtx<float> mX = "[ 8 6 2; 7 4 9; 4 8 2 ]";
    mX.inverse().print();

    // Make a random matrix
    Mtx<double> mA(20,5);
    Rand rnd;
    mA.random_gaussian(rnd);

    // Compute the SVD
    MatrixDecompSVD<double> svd(mA);

    // Compute the reconstruction error
    double maxerror = (svd.matrixU() * svd.vectorS().diag_matrix() * svd.matrixV().t() - mA).max_abs();

    // Fill upper left 4x4 region of A with zero
    mA.block(0,0,4,4) = 0.0;

    // Compute the dot product of two columns of A
    double dotprod = mA.col(2).dot_product(mA.col(3));
    
    // Copy the first 5 elements of the first column vector in U over to the diagonal of a sub
    // block of A.
    // Its important to understand that these operations manipulate pointers to views of memory
    // and the copy_from() function is the only one that actually does any data copying.
    mA.block(10,0,5,5).diag().copy_from(svd.matrixU.col(0).head(5));

### Build Instructions
This code was developed on the Mac. We will be adding some make files eventually, but the code itself is not very demanding to build. To compile against it, merely include the imp_core.h file. All the functions and objects are within the 'im' namespace.

