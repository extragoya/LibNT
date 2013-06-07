cmake_minimum_required (VERSION 2.6)

set (ASSIGN_SOURCES mia_assign_test.cpp)
add_executable(MIAAssignTests ${ASSIGN_SOURCES})
ADD_TEST(MIAAssignTests ${CMAKE_CURRENT_BINARY_DIR}/MIAAssignTests)

set (CONSTR_SOURCES mia_constr_test.cpp)
add_executable(MIAConstrTests ${CONSTR_SOURCES})
ADD_TEST(MIAConstrTests ${CMAKE_CURRENT_BINARY_DIR}/MIAConstrTests)

set (MERGE_SOURCES mixed_mia_merge_test.cpp)
add_executable(MixedMIAMergeTests ${MERGE_SOURCES})
ADD_TEST(MixedMIAMergeTests ${CMAKE_CURRENT_BINARY_DIR}/MixedMIAMergeTests)

set (MULT_SOURCES mixed_mia_mult_test.cpp)
add_executable(MixedMIAMultTests ${MULT_SOURCES})
ADD_TEST(MixedMIAMultTests ${CMAKE_CURRENT_BINARY_DIR}/MixedMIAMultTests)


if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(MIAAssignTests ${Boost_LIBRARIES})
	target_link_libraries(MIAConstrTests ${Boost_LIBRARIES})
	target_link_libraries(MixedMIAMergeTests ${Boost_LIBRARIES})
	target_link_libraries(MixedMIAMultTests ${Boost_LIBRARIES})
endif()