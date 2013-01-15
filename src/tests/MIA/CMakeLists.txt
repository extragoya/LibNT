cmake_minimum_required (VERSION 2.6)

set (ASSIGN_SOURCES mia_assign_test.cpp)
add_executable(MIAAssignTests ${ASSIGN_SOURCES})
ADD_TEST(MIAAssignTests ${CMAKE_CURRENT_BINARY_DIR}/MIAAssignTests)


