# ImpressiveAI
Impressive Machines's C++ Library for Machine Learning and Computer Vision

The aim of this library is to create a self-contained set of tools for many applications in machine learning and computer vision. It is broken into a number of sections which are useful in isolation. The main part of the library is "core" which is a library for linear algebra and numberical computation. This part currently has a lot of functionality. Other parts are under development and will be added in the future.

## Core
Core provides many common tools for linear algebra and is built on very flexible matrix and vector objects. The features include:

* Flexible and light-weight matrix and vector views which provide a way to refer to one and two dimensional memory. They have arbitrary row and columns strides and so can operate on many data formats and with any data type.
* Optimized low-level matrix block functions which work on vector and matrix types to implement common operations, such as math operations, data sorting, conversion, statistics and reduce operations.
* A full suite of conforming BLAS operations supporting BLAS 1, 2 and 3 linear algebra manipulations. These make use of matrix and vector views for their arguments.
* Matrix and vector types which manage shared memory and can be used to express quite elaborate mathematical operations neatly and efficiently.
* Real and complex types are supported for many math operations.
* Typical matrix decompositions:
    * LU decomposition.
    * Cholesky LDLT and LLT decomposition.
    * QR decomposition with and without pivoting.
    * Tridiagonal decomposition.
    * Eigenvalue decomposition for symmetric matrices.
    * Jacobi SVD with QR preconditioning.
    * Special case 2x2 eigen and SVD decompositions.
* Efficient convolution for 1D and 2D signals.
* Matrices and vectors can be saved and printed in a variety of formats, including matlab text format.
* Convenient interfaces for random number generation.
* Solution to quadratic and cubic polynomials.
* Principal components analysis.
* Estimation of covariances.
* Generation of multidimensioanl Gaussian noise with specified covariance.
* Computation of multivariate Gaussians and Mahalanobis distance.
* Arbitrary radix FFT in 1D and 2D (TBD).
* Minimization of functions (TBD).
* Regression (TBD).
* Tensors (TBD).

### Dependencies
There are no external dependencies for the core library.

### Philosophy
There are many linear algebra libraries available that are built using C++. Many of the popular ones make extensive use of template metaprogramming to increase efficiency. This library takes a different approach. Not many developers really understand C++ templates well and so complex template code is hard to debug and maintain. This library aims towards simplicity of code design. It does this without giving up too much run time efficiency by making extensive use of light weight references to matrix memory. It avoids un-neccesary data copies by using shared pointers where possible. The main goal was to make the code easy to understand and still quite fast. The primary source for many of the algorithms was the GNU scientific library and so the routines are generally expected to be of good quality and have performed well in testing so far.

### Examples
    // Make a random matrix and compute the SVD and then compute the reconstruction error
    im::Mtx<double> mA(20,5);
    im::Rand rnd;
    mA.random_gaussian(rnd);
    im::MatrixDecompSVD<double> svd(mA);
    double maxerror = (svd.matrixU() * svd.vectorS().diag_matrix() * svd.matrixV().t() - m2).max_abs();

    // Fill upper left 4x4 region of A with zero
    mA.block(0,0,4,4) = 0.0;

    // Compute the dot product of two columns of A
    double dotprod = mA.col(2).dot_product(mA.col(3));
    
    // copy the first 5 elements of the first column vector in U over to the diagonal of a sub block of A.
    // Its important to understand that these operations manipulate pointers to views of memory and the
    // copy_from() function is the only one that actually does any data copying.
    mA.block(10,0,5,5).diag().copy_from(svd.matrixU.col(0).head(5));


