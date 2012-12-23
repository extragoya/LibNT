cmake_minimum_required (VERSION 2.6)

project(SparseLattice_example)
set (SOURCES slattice_1.cpp slattice_test.h)

add_executable(SparseLattice_example ${SOURCES})

INCLUDE_DIRECTORIES(${SuperLU_INCLUDE})
TARGET_LINK_LIBRARIES(SparseLattice_example ${blas_LIB}
                      ${SuperLU_LIB}
)


