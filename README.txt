LibNT is header only (except for MATLAB binaries), so using it should be fairly simple.

************Note on name change************
The library is currently transferring from its old name, LibMIA, to its new name LibNT. Thus, you will likely encounter nomenclature referring to multi-index arrays (MIAs) instead of numeric tensors (NTs). Eventually all references and class names using MIA will be removed. 


************LibNT: C************

Requirements to use LibNT within your own C++ program. You don't need to build anything, you just need to include headers from two other header-only libraries.

1. Eigen http://eigen.tuxfamily.org/: latest version LibNT has been tested on is Eigen 3.2.5. Older versions will not work, as Eigen's previous sparseQR and sparseCholesky didn't work with outside data. Include the root directory within your build.
    -Eigen has optional support for the Intel MKL library, enabling faster execution of dense operations, as explained here: http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html. If have the libraries and you want to enable this you must define the appropriate macro (see link) and also link to the appropriate MKL libraries.
2. Boost http://www.boost.org/: you only need the headers, don't bother compiling any of the binary libraries, unless you have a need to/want to. Include the Boost root directory within your include directories.
    -An exception is if you want to build LibNT's tests. In this case you can optionally link to boost's unit_test_framework library to speed up compilation. The easiest way to build the tests is to use the CMake build system configured using the root CMakeLists.txt in LibNT's source directory. You can follow the MATLAB instructions below to do this.
3. Include the LibMIA headers found in <libmia_root>/src/multi_index_array.
4. A C++11 compliant compiler. As of writing MSVC 2015 has just been released and should meet the requirements. Intel and GCC definitely work. 

Note: an example file can be found in <libmia_root>/example with a CMake configuration file. Some example code can be found there. In addition the CMakeLists.txt file there details some of the build options you can use.

************NTToolbox: MATLAB Tensor classes with an Interface to LibNT's Algorithms************

--------
First, like above, you need to obtain the following two header-only libraries and a C++11 compliant compiler
1. Eigen http://eigen.tuxfamily.org/: latest version LibNT has been tested on is Eigen 3.2.5. Older versions will not work, as Eigen's previous sparseQR and sparseCholesky didn't work with outside data. Include the root directory within your build.   
2. Boost http://www.boost.org/: you only need the headers, don't bother compiling any of the binary libraries, unless you have a need to/want to. Include the Boost root directory within your include directories.
3. A C++11 compliant compiler. As of writing MSVC 2015 has just been released and should meet the requirements. Intel and GCC definitely work. 
--------    

The recommended way to build the MATLAB mex files is through the CMake build system http://www.cmake.org/. If using this option, download the cmake build system. Point the root directory of LibNT as the source directory. Choose a build location. Configure the build for your preferred build system, e.g., gcc or MSFT Visual Studio. Then generate the build. You will need to set the following options.

1. Eigen: point Eigen_ROOT to the root directory of the Eigen library
2. Boost: point Boost_ROOT to the root directory of the Boost library. 
3. <Optional> If you want to perform the C++ tests, turn MIA_BUILD_TESTS on.
   - If you've built Boost's unit_test_framework library, you can turn MIA_USE_HEADER_ONLY_TESTS off to link to that library, speeding up the compilation of the tests. You may need to also point Boost_LIBRARYDIR to the location of your compiled boost libraries
4. Turn MIA_BUILD_MATLAB on to build the MEX interface. 
5. Point MATLAB_INCLUDE to the MATLAB header directory, usually <MATLAB_install>/extern/include. Point MATLAB_LIB_DIR to the MATLAB library directory, e.g., <MATLAB_install>extern/lib/win64/microsoft
6. <Optional> If you have Intel MKL, you can use the library to speed up execution of Eigen's routines. Turn USE_MKL on and point MKL_ROOT to the root directory of the Intel MKL installation.
5. You should be able to now successfully generate your build using CMake. Open up the generated project (or if using command line, navigate to the build directory). Build the MATLAB project, this will build your .mex files. 
6. Add the MATLAB classes in <libnt_root>/src/MATLAB to MATLAB's search path. Also add your build directory for the .mex files to MATLAB's search path, i.e., <libnt_BUILD_root>/src/MATLAB/mex. Now you should be good to go to use the NT classes. 


Instead of using the CMake build system, these options can all be manually performed, but you will have to perform several extra steps. You can consult the CMakeLists.txt files in the LibNT source for guidance on this. But, the CMake build system remains the best option. 
