# Metaphor
C++ Library for Machine Learning and Computer Vision
Impressive Machines LLC

The eventual aim of this library is to create a self-contained set of tools for many applications in machine learning and computer vision. It is broken into a number of sections which are useful in isolation. The main part of the library is "core" which is a library for linear algebra and numerical computation. This part currently has a lot of functionality. Other parts are under development and will be added in the future.

## Core
Core provides many common tools for linear algebra and statistics and is built on very flexible matrix and vector objects. The features include:

* General purpose light-weight matrix and vector view classes:
    * View objects can reference existing data with arbitrary signed row and column stride.
    * Any data type supported.
    * Many common data operations, such as sub-views copying, re-ordering, and decimation.
    * Matrices and vectors can be saved and printed in a variety of formats, including Matlab text format.
* Linear algebra operations for float, double, and complex types:
    * A full suite of optimized standard BLAS operations supporting BLAS 1, 2 and 3.
    * Optimized low-level matrix block math functions which work on vector and matrix view types.
    * Math operations include data sorting, conversion, statistics, and reduce operations.
* Matrix and vector memory objects
    * Tightly integrated with view objects
    * Offer a multitude of member functions for common linear algebra operations.
    * Use shared pointers to memory to avoid data copying.
* Typical matrix decompositions:
    * LU decomposition.
    * Cholesky LDLT and LLT decomposition.
    * QR decomposition with and without pivoting.
    * Tridiagonal decomposition.
    * Eigenvalue decomposition for symmetric matrices.
    * Jacobi SVD with QR preconditioning.
    * Special case 2x2 Eigen and SVD decompositions.
* Signal processing:
    * Fast block-based convolution for 1D and 2D signals.
    * Arbitrary radix FFT in 1D and 2D.
* Statistics:
    * Convenient interfaces for random number generation.
    * Principal components analysis.
    * Estimation of covariances.
    * Generation of multidimensional Gaussian noise with specified covariance.
    * Computation of multivariate Gaussians and Mahalanobis distance.
* Solution to quadratic and cubic polynomials.
* Quadratic fitting and interpolation in 1D and 2D.
* Least squares and robust 2D line fitting.
* Multi-dimensional regression.
* Geometric 2D and 3D rotations and transforms.
* Quaternions with SLERP and rotation matrix / euler / axis-angle integration.
* Numerical root finding for scalar functions.
* Optimization of functions:
    * Bracketing, line minimization, and line minimization with derivatives.
    * Levenberg Marquardt minimization.
    * Stochastic gradients with optional momentum or acceleration.
    * Powell minimization.
    * Conjugate gradients.
    * BFGS.
    * Limited memory BFGS.

The matrix and vector view objects have independent signed row and columns strides and so they can operate on many data formats and with any data type - they do not require classical row or column major memory organization. All operations support positive or negative row or column strides. Therefore they can refer to memory which is reversed, decimated, or fragmented, such as all red pixels of an RGB image with arbitrary row stride such as being stored in bottom up order. 

Decimating, reversing, or extracting rows, columns, blocks, or diagonals of matrices is therefore a constant time operation. Most operations in the library involve passing only the small view structures as arguments to functions (often by reference) and unnecessary data copying is avoided. Matrix and vector classes exist to dynamically create memory blocks where needed and these make use of standard shared pointers. Alternatively, memory can be statically allocated by the user and trivially wrapped into matrix view objects for use with any library functions. Other data blocks in related libraries, such as images, may be referenced easily thought this library's matrix or vector view classes.

### Dependencies
There are no external dependencies for the core library.

### Philosophy
There are many linear algebra libraries available that are built using C++. Many of the popular ones make extensive use of template meta-programming to increase efficiency. This library takes a different approach. Not many developers really understand C++ templates well and so complex template code is hard to debug and maintain. 

This library aims towards simplicity of code design. It does this without giving up too much run time efficiency by making extensive use of light weight references to matrix memory. It avoids unneccesary data copies by using shared pointers where possible. The main goal was to make the code easy to understand and still quite fast.

### Examples
    #include "meta_core.h"
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
1. Download and install CMAKE - you will need version 3.1 or above
2. Create a build directory wherever you want to build the library
3. From the command line, cd to your build directory
4. Run cmake "path/to/the/metaphor/git/source/directory"
5. If all is well, run make, which should compile everything and generate a library file
6. Write your code and refer to metaphor objects using namespace "im"
6. Include the file meta_core.h for the core library at the start of your code





