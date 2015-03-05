cmake_minimum_required (VERSION 2.6)
PROJECT(LibMIA_example)
#turn on C++11 flag
if(CMAKE_COMPILER_IS_GNUCXX OR MINGW OR CYGWIN)
	#gcc accepts a tag for restricted pointers (https://gcc.gnu.org/onlinedocs/gcc/Restricted-Pointers.html), other compilers have their own flags and rules, so you can adjust accordingly.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x ${CMAKE_CXX_FLAGS} -Drestrict_libmia=__restrict__")  
  
  #You can also add your own more aggressive optimizations here if you wish.
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline-functions -ftemplate-depth-128") 
endif()

#point CMake to your boost directory
SET(BOOST_ROOT CACHE PATH /path/to/boost)
find_package(Boost REQUIRED)#EXACT COMPONENTS unit_test_framework timer system chrono)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS}) 
endif()

#point CMake to your Eigen directory
SET(Eigen_ROOT  CACHE PATH /path/to/eigen)
include_directories(${Eigen_ROOT}) 
#point CMake to your <libmia_root>/src/multi_index_array  directory
SET(LibMIA_ROOT  CACHE PATH /path/to/libmia_root/src/multi_index_array)
include_directories(${LibMIA_ROOT}) 


#create a target executable called example_cmake using the main.cpp source file 
add_executable(example_cmake main.cpp)
