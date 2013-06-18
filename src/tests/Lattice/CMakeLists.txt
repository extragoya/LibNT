cmake_minimum_required (VERSION 2.8)




set (CONSTR_SOURCES latticeconstrtests.cpp)
add_executable(LatticeConstrTests ${CONSTR_SOURCES})
ADD_TEST(LatticeConstrTests ${CMAKE_CURRENT_BINARY_DIR}/LatticeConstrTests)

set (MIXED_MULT_SOURCES mixedlatticemulttests.cpp)
add_executable(MixedLatticeMultTests ${MIXED_MULT_SOURCES})
ADD_TEST(MixedLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/MixedLatticeMultTests)



if(LIBMIA_USE_SPARSE_SOLVE)
	set (MIXED_SOLVE_SOURCES mixed_lattice_solve_test.cpp)
	add_executable(MixedLatticeSolveTests ${MIXED_SOLVE_SOURCES})
	ADD_TEST(MixedLatticeSolveTests ${CMAKE_CURRENT_BINARY_DIR}/MixedLatticeSolveTests)
	if(NOT MIA_USE_HEADER_ONLY_TESTS)
		target_link_libraries(MixedLatticeSolveTests ${SuperLU_LIB} ${blas_LIB} ${Boost_LIBRARIES})
	else()
		target_link_libraries(MixedLatticeSolveTests ${SuperLU_LIB} ${blas_LIB})
	endif()
endif(LIBMIA_USE_SPARSE_SOLVE)

if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(MixedLatticeMultTests ${Boost_LIBRARIES})
	target_link_libraries(LatticeConstrTests ${Boost_LIBRARIES})
endif()