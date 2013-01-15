cmake_minimum_required (VERSION 2.6)


SET(Boost_ROOT  CACHE PATH /path/to/boost)
SET(Eigen_ROOT  CACHE PATH /path/to/eigen)
SET(USE_SPARSE_SOLVE TRUE CACHE BOOL "build routines requiring sparse solver libs")

SET(blas_LIB  CACHE FILEPATH /path/to/blas_lib)
SET(SuperLU_LIB  CACHE FILEPATH /path/to/superlu_lib)
SET(SuperLU_INCLUDE  CACHE PATH /path/to/superlu_headers)

INCLUDE_DIRECTORIES(${Boost_ROOT})
INCLUDE_DIRECTORIES(${Eigen_ROOT})
INCLUDE_DIRECTORIES(${SuperLU_INCLUDE})

add_subdirectory(multi_index_array)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/multi_index_array)
add_subdirectory(tests)
