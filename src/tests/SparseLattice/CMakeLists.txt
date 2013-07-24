cmake_minimum_required (VERSION 2.8)

#set (COMPARE_SOURCES sparse_mia_compare_test.cpp)
#add_executable(SparseMIACompareTests ${COMPARE_SOURCES})
#set (FUNCTION_SOURCES sparse_mia_functions.cpp)
#add_executable(SparseMIAFunctionTests ${FUNCTION_SOURCES})

set (MULT_SOURCES sparse_lattice_mult_test.cpp)
set (SOLVE_SOURCES sparse_lattice_solve_test.cpp)
add_executable(SparseLatticeMultTests ${MULT_SOURCES})
add_executable(SparseLatticeSolveTests ${SOLVE_SOURCES})
if(NOT MIA_USE_HEADER_ONLY_TESTS)
	 #target_link_libraries(SparseMIACompareTests ${Boost_LIBRARIES})
	 #target_link_libraries(SparseMIAFunctionTests ${Boost_LIBRARIES})	 
	 target_link_libraries(SparseLatticeMultTests ${Boost_LIBRARIES})
	 target_link_libraries(SparseLatticeSolveTests ${Boost_LIBRARIES})
endif()
ADD_TEST(SparseLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeMultTests)
ADD_TEST(SparseLatticeSolveTests ${CMAKE_CURRENT_BINARY_DIR}/SparseLatticeSolveTests)


#ADD_TEST(SparseMIACompareTests ${CMAKE_CURRENT_BINARY_DIR}/SparseMIACompareTests)
#ADD_TEST(SparseMIAFunctionTests ${CMAKE_CURRENT_BINARY_DIR}/SparseMIAFunctionTests)




