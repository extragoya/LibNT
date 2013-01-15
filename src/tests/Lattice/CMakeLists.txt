cmake_minimum_required (VERSION 2.6)


set (MULT_SOURCES latticemulttests.cpp)
add_executable(DenseLatticeMultTests ${MULT_SOURCES})
ADD_TEST(DenseLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/DenseLatticeMultTests)

set (CONSTR_SOURCES latticeconstrtests.cpp)
add_executable(LatticeConstrTests ${CONSTR_SOURCES})
ADD_TEST(LatticeConstrTests ${CMAKE_CURRENT_BINARY_DIR}/LatticeConstrTests)


file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

