cmake_minimum_required (VERSION 2.6)

set (ASSIGN_SOURCES mia_assign_test.cpp)
add_executable(MIAAssignTests ${ASSIGN_SOURCES})
ADD_TEST(MIAAssignTests ${CMAKE_CURRENT_BINARY_DIR}/MIAAssignTests)

set (CONSTR_SOURCES mia_constr_test.cpp)
add_executable(MIAConstrTests ${CONSTR_SOURCES})
ADD_TEST(MIAConstrTests ${CMAKE_CURRENT_BINARY_DIR}/MIAConstrTests)
if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(MIAAssignTests ${Boost_LIBRARIES})
	target_link_libraries(MIAConstrTests ${Boost_LIBRARIES})
endif()