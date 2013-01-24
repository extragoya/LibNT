cmake_minimum_required (VERSION 2.6)



SET(Eigen_ROOT  CACHE PATH /path/to/eigen)
#if you want to use sparse solution of equations, you need to have built Super_LU, and then you must point cmake to the blas_LIB, SuperLU_LIB and Super_LU include directories
SET(USE_SPARSE_SOLVE TRUE CACHE BOOL "build routines requiring sparse solver libs")
#only change this flag if you've built the boost::test libs
SET(MIA_USE_HEADER_ONLY_TESTS TRUE CACHE BOOL "use the single header option in boost::test, which will work universally but will compile slowly")

SET(BOOST_ROOT  CACHE PATH /path/to/blas_lib)

SET(blas_LIB  CACHE FILEPATH /path/to/blas_lib)
SET(SuperLU_LIB  CACHE FILEPATH /path/to/superlu_lib)
SET(SuperLU_INCLUDE  CACHE PATH /path/to/superlu_headers)

#if we're using the single-header option in boost::test, only include the boost directories
if(MIA_USE_HEADER_ONLY_TESTS)
	
	find_package(Boost 1.52.0 REQUIRED)
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
		
	endif()
else(MIA_USE_HEADER_ONLY_TESTS) #otherwise, try to find the precompiled libs for the unit_test_framework
	
	set(Boost_USE_STATIC_LIBS        ON)
	set(Boost_USE_MULTITHREADED      ON)
	set(Boost_USE_STATIC_RUNTIME    OFF)
	find_package(Boost 1.52.0 REQUIRED COMPONENTS unit_test_framework )
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
	endif()
endif(MIA_USE_HEADER_ONLY_TESTS)

#setup the config file for LibMIA
configure_file (
  "MIAConfig.h.in"
  "${PROJECT_BINARY_DIR}/MIAConfig.h"
  )
include_directories("${PROJECT_BINARY_DIR}")



INCLUDE_DIRECTORIES(${Eigen_ROOT})
if(USE_SPARSE_SOLVE)
	INCLUDE_DIRECTORIES(${SuperLU_INCLUDE})
endif()

add_subdirectory(multi_index_array)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/multi_index_array)
add_subdirectory(tests)
