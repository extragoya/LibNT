cmake_minimum_required (VERSION 2.6)

set (MULT_SOURCES dense_mia_mult_test.cpp)
add_executable(DenseMIAMultTests ${MULT_SOURCES})
ADD_TEST(DenseMIAMultTests ${CMAKE_CURRENT_BINARY_DIR}/DenseMIAMultTests)


