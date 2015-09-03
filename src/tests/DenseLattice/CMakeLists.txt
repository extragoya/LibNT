cmake_minimum_required (VERSION 2.6)


set (MULT_SOURCES denselatticemulttests.cpp)
add_executable(DenseLatticeMultTests ${MULT_SOURCES})
ADD_TEST(DenseLatticeMultTests ${CMAKE_CURRENT_BINARY_DIR}/DenseLatticeMultTests)



file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(DenseLatticeMultTests ${Boost_LIBRARIES})
endif()