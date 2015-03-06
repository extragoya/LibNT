LibMIA is header only (except for MATLAB binaries), so using it should be fairly simple.

**C++ Interface**

Requirements to use LibMIA within your own C++ program. You don't need to build anything, you just need to include headers from two other header-only libraries.

1. Eigen: latest version it has been tested on is Eigen 3.2.1. Then, please follow the instructions in the README file of <libmia_root>/Eigen. This directory contains modified Eigen source files that are needed for the execution of LibMIA. Include the Eigen root directory within your include directories.
2. Boost: you should only need the headers, don't bother compiling any of the binary libraries, unless you have a need to/want to. Include the Boost root directory within your include directories.
3. Include the LibMIA headers in <libmia_root>/src/multi_index_array.
4. A C++11 compliant compiler. As of writing MSVC still doesn't meet the requirements, but its upcoming 2015 release should.

Note: an example file can be found in <libmia_root>/example with a CMake configuration file.

**MATLAB Interface**

In this case Requirements to build MATLAB mex files (need CMake).

1. Eigen: latest version it has been tested on is Eigen 3.2.1. Then, please follow the instructions in the README file of <libmia_root>/Eigen. This directory contains modified Eigen source files that are needed for the execution of LibMIA. 
2. Boost: you should only need the headers, don't bother compiling any of the binary libraries, unless you have a need to/want to.
3. MATLAB installation. Include the MATLAB header directory within your include directories, usually <MATLAB_install>/extern/include. Include the MATLAB library directory within your library search directories, usually <MATLAB_install>extern/lib/win64/microsoft
4. If using CMake (recommended), the LibMIA CMake configuration includes options for 1, 2, and 3. Also make sure you turn on the MIA_BUILD_MATLAB option.
5. Build the MATLAB project, this will build your .mex files. 
6. Add the MATLAB classes in <libmia_root>/src/MATLAB to MATLAB's search path. Also add your build directory for the .mex files to MATLAB's search path. Now you should be good to go to use the MIA classes. 




