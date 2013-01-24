cmake_minimum_required (VERSION 2.6)
list( APPEND CMAKE_CXX_FLAGS "-std=gnu++0x ${CMAKE_CXX_FLAGS}")
PROJECT(LibMIA)
set (LibMIA_VERSION_MAJOR 0)
set (LibMIA_VERSION_MINOR 2)
add_subdirectory (src)

