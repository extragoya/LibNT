cmake_minimum_required (VERSION 2.6)


set (SOURCES latticemulttests.cpp)
add_executable(DenseLatticeMultTests ${SOURCES})
ADD_TEST(DenseLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/DenseLatticeMultTests)


file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

