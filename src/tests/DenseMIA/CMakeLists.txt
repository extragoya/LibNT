cmake_minimum_required (VERSION 2.6)

set (MULT_SOURCES dense_mia_mult_test.cpp)
add_executable(DenseMIAMultTests ${MULT_SOURCES})
if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(DenseMIAMultTests ${Boost_LIBRARIES})
endif()
ADD_TEST(DenseMIAMultTests ${CMAKE_CURRENT_BINARY_DIR}/DenseMIAMultTests)


