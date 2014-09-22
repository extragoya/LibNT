cmake_minimum_required (VERSION 2.8)



SET(Eigen_ROOT  CACHE PATH /path/to/eigen)

#only change this flag if you've built the boost::test libs
SET(MIA_USE_HEADER_ONLY_TESTS TRUE CACHE BOOL "use the single header option in boost::test, which will work universally but will compile slowly")

SET(MIA_BUILD_TESTS FALSE CACHE BOOL "if you would like to build tests")

SET(MIA_BUILD_MATLAB TRUE CACHE BOOL "if you would like to build MATLAB mex files")

SET(BOOST_ROOT  CACHE PATH /path/to/boost)
set(Boost_DIR ${BOOST_ROOT})





if(NOT MIA_USE_HEADER_ONLY_TESTS) #otherwise, try to find the precompiled libs for the unit_test_framework
	set(Boost_USE_STATIC_LIBS        ON)
	set(Boost_USE_MULTITHREADED      OFF)
	set(Boost_USE_STATIC_RUNTIME    OFF)
	set (Boost_USE_DEBUG_RUNTIME OFF)
	find_package(Boost REQUIRED COMPONENTS unit_test_framework timer system chrono)
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
		message(${Boost_LIBRARIES})
	endif()
	
else(NOT MIA_USE_HEADER_ONLY_TESTS) #if we're using the single-header option in boost::test, only include the boost directories
	find_package(Boost REQUIRED)
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
		
	endif()
	
endif(NOT MIA_USE_HEADER_ONLY_TESTS)
set(Boost_install_include "${Boost_INCLUDE_DIRS}" PARENT_SCOPE)

#setup the config file for LibMIA
configure_file (
  "MIAConfig.h.in"
  "${PROJECT_BINARY_DIR}/MIAConfig.h"
  )
include_directories("${PROJECT_BINARY_DIR}")



INCLUDE_DIRECTORIES(${Eigen_ROOT})


add_subdirectory(multi_index_array)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/multi_index_array)
if(MIA_BUILD_TESTS)
	add_subdirectory(tests)
endif(MIA_BUILD_TESTS)

if(MIA_BUILD_MATLAB)
	add_subdirectory(MATLAB)
endif(MIA_BUILD_MATLAB)
