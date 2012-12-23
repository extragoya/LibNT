LibMIA 0.1

Very Early Release:
Dense and Sparse Lattice support
Early MIA expression checking and dispatching

Dependencies:

C++11 Compiler that supports __COUNTER__ macro - ie gcc 4.3 and later and MSVC

Almost all dependencies are header-based
Download Boost http://www.boost.org/ and included it within your search directories
Download Eigen 3 http://eigen.tuxfamily.org and include it within your search directories

The exception is if you're solving sparse equations
Sparse solution of equations requires the SuperLU extension to Eigen. This requires generating the SuperLU library and linking to it. Compiling the SuperLU library requires in turn a blas library.
Note that this also relies on a patch APH made to Eigen. A patch has been suggested, the status of it can be tracked http://eigen.tuxfamily.org/bz/show_bug.cgi?id=454


