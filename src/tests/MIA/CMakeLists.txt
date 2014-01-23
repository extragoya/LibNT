cmake_minimum_required (VERSION 2.6)

set (ASSIGN_SOURCES mia_assign_test.cpp)
add_executable(MIAAssignTests ${ASSIGN_SOURCES})
ADD_TEST(MIAAssignTests ${CMAKE_CURRENT_BINARY_DIR}/MIAAssignTests)

set (MIXED_ASSIGN_SOURCES mixed_mia_assign_test.cpp)
add_executable(MixedMIAAssignTests ${MIXED_ASSIGN_SOURCES})
ADD_TEST(MixedMIAAssignTests ${CMAKE_CURRENT_BINARY_DIR}/MixedMIAAssignTests)

set (UNARY_SOURCES mia_unary_test.cpp)
add_executable(MIAUnaryTests ${UNARY_SOURCES})
ADD_TEST(MIAUnaryTests ${CMAKE_CURRENT_BINARY_DIR}/MIAUnaryTests)

set (CONSTR_SOURCES mia_constr_test.cpp)
add_executable(MIAConstrTests ${CONSTR_SOURCES})
ADD_TEST(MIAConstrTests ${CMAKE_CURRENT_BINARY_DIR}/MIAConstrTests)

set (MERGE_SOURCES mixed_mia_merge_test.cpp)
add_executable(MixedMIAMergeTests ${MERGE_SOURCES})
ADD_TEST(MixedMIAMergeTests ${CMAKE_CURRENT_BINARY_DIR}/MixedMIAMergeTests)

set (MULT_SOURCES mixed_mia_mult_test.cpp)
add_executable(MixedMIAMultTests ${MULT_SOURCES})
ADD_TEST(MixedMIAMultTests ${CMAKE_CURRENT_BINARY_DIR}/MixedMIAMultTests)

#need to do separate tests for float and double because the of section limits in obj files when compiling using MinGw. 
set (SOLVE_SOURCES_FLOAT mixed_mia_solve_test.hpp mixed_mia_solve_test_float.cpp)
add_executable(MixedMIASolveTestsFloat ${SOLVE_SOURCES_FLOAT})
ADD_TEST(MixedMIASolveTestsFloat ${CMAKE_CURRENT_BINARY_DIR}/MixedMIASolveTestsFloat)

set (SOLVE_SOURCES_DOUBLE mixed_mia_solve_test.hpp mixed_mia_solve_test_double.cpp)
add_executable(MixedMIASolveTestsDouble ${SOLVE_SOURCES_DOUBLE})
ADD_TEST(MixedMIASolveTestsDouble ${CMAKE_CURRENT_BINARY_DIR}/MixedMIASolveTestsDouble)

if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(MIAAssignTests ${Boost_LIBRARIES})
	target_link_libraries(MixedMIAAssignTests ${Boost_LIBRARIES})
	target_link_libraries(MIAUnaryTests ${Boost_LIBRARIES})
	target_link_libraries(MIAConstrTests ${Boost_LIBRARIES})
	target_link_libraries(MixedMIAMergeTests ${Boost_LIBRARIES})
	target_link_libraries(MixedMIAMultTests ${Boost_LIBRARIES})
	target_link_libraries(MixedMIASolveTestsFloat ${Boost_LIBRARIES})
	target_link_libraries(MixedMIASolveTestsDouble ${Boost_LIBRARIES})
endif()

