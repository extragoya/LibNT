cmake_minimum_required (VERSION 2.6)
PROJECT(LibMIA_example)
#turn on C++11 flag
list( APPEND CMAKE_CXX_FLAGS "-std=gnu++0x ${CMAKE_CXX_FLAGS}")
#find the LibMIA package (assumes you have used CMake to create a build and installation, otherwise find_package will not work)
find_package(LibMIA REQUIRED)
#include the LibMIA headers in the project
include_directories(${LIBMIA_INCLUDE_DIRS}) 
#create a target executable called example_cmake using the main.cpp source file 
add_executable(example_cmake main.cpp)
