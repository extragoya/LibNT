cmake_minimum_required (VERSION 2.8)



set (MULT_SOURCES sparse_lattice_mult_test.cpp)
set (MULT_SOURCES_SPECIAL sparse_lattice_mult_test_specialized.cpp)
set (SOLVE_SOURCES_FLOAT sparse_lattice_solve_test.hpp sparse_lattice_solve_test_float.cpp)
set (SOLVE_SOURCES_DOUBLE sparse_lattice_solve_test.hpp sparse_lattice_solve_test_double.cpp)
SET (FUNCTION_SOURCES sparse_lattice_functions.cpp)
add_executable(SparseLatticeMultTests ${MULT_SOURCES})
add_executable(SparseLatticeMultTestsSpecial ${MULT_SOURCES_SPECIAL})
add_executable(SparseLatticeSolveTestsFloat ${SOLVE_SOURCES_FLOAT})
add_executable(SparseLatticeSolveTestsDouble ${SOLVE_SOURCES_DOUBLE})
add_executable(SparseLatticeFunctions ${FUNCTION_SOURCES})
if(NOT MIA_USE_HEADER_ONLY_TESTS)	 
	 target_link_libraries(SparseLatticeMultTests ${Boost_LIBRARIES})
	 target_link_libraries(SparseLatticeMultTestsSpecial ${Boost_LIBRARIES})
	 target_link_libraries(SparseLatticeSolveTestsFloat ${Boost_LIBRARIES})
	 target_link_libraries(SparseLatticeSolveTestsDouble ${Boost_LIBRARIES})
	 target_link_libraries(SparseLatticeFunctions ${Boost_LIBRARIES})
endif()
ADD_TEST(SparseLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeMultTests)
ADD_TEST(SparseLatticeMultTestsSpecial ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeMultTestsSpecial)
ADD_TEST(SparseLatticeSolveTestsFloat ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeSolveTestsFloat)
ADD_TEST(SparseLatticeSolveTestsDouble ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeSolveTestsDouble)
ADD_TEST(SparseLatticeFunctions ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeFunctions)






