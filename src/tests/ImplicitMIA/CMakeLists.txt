cmake_minimum_required (VERSION 2.6)


set (BLOCK_SOURCES implicit_mia_block_functions.cpp)
add_executable(ImplicitMIABlockTests ${BLOCK_SOURCES})
set (FUNCTION_SOURCES implicit_mia_functions.cpp)
add_executable(ImplicitMIAFunctionTests ${FUNCTION_SOURCES})
if(NOT MIA_USE_HEADER_ONLY_TESTS)
	target_link_libraries(ImplicitMIAFunctionTests ${Boost_LIBRARIES})
	target_link_libraries(ImplicitMIABlockTests ${Boost_LIBRARIES})	
endif()
ADD_TEST(ImplicitMIAFunctionTests ${CMAKE_CURRENT_BINARY_DIR}/ImplicitMIAFunctionTests)
ADD_TEST(ImplicitMIABlockTests ${CMAKE_CURRENT_BINARY_DIR}/ImplicitMIABlockTests)



