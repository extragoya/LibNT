cmake_minimum_required (VERSION 2.8)



SET(Eigen_ROOT  CACHE PATH /path/to/eigen)

#only change this flag if you've built the boost::test libs
SET(MIA_USE_HEADER_ONLY_TESTS TRUE CACHE BOOL "use the single header option in boost::test, which will work universally but will compile slowly")

SET(BOOST_ROOT  CACHE PATH /path/to/boost)
set(Boost_DIR ${BOOST_ROOT})




#if we're using the single-header option in boost::test, only include the boost directories
if(MIA_USE_HEADER_ONLY_TESTS)
	
	find_package(Boost 1.53.0 REQUIRED EXACT)
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
		
	endif()
else(MIA_USE_HEADER_ONLY_TESTS) #otherwise, try to find the precompiled libs for the unit_test_framework
	
	set(Boost_USE_STATIC_LIBS        ON)
	set(Boost_USE_MULTITHREADED      ON)
	set(Boost_USE_STATIC_RUNTIME    OFF)
	find_package(Boost 1.53.0 REQUIRED EXACT COMPONENTS unit_test_framework timer system chrono )
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
	endif()
endif(MIA_USE_HEADER_ONLY_TESTS)
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
add_subdirectory(tests)
