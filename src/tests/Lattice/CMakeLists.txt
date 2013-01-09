cmake_minimum_required (VERSION 2.6)


set (MULT_SOURCES latticemulttests.cpp)
add_executable(DenseLatticeMultTests ${MULT_SOURCES})
ADD_TEST(DenseLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/DenseLatticeMultTests)

set (CONSTR_SOURCES latticeconstrtests.cpp)
add_executable(DenseLatticeConstrTests ${CONSTR_SOURCES})
ADD_TEST(DenseLatticeConstrTests ${CMAKE_CURRENT_BINARY_DIR}/DenseLatticeConstrTests)


file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

