cmake_minimum_required (VERSION 2.8)

SET(MATLAB_INCLUDE  CACHE PATH /path/to/matlab_include)
SET(MATLAB_LIB_DIR  CACHE PATH /path/to/matlab_lib_dir)
include_directories(${MATLAB_INCLUDE}) 


add_subdirectory(mex)
