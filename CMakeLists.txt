cmake_minimum_required (VERSION 2.6)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
macro (add_sources)
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND SRCS "${_relPath}/${_src}")
        else()
            list (APPEND SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (SRCS ${SRCS} PARENT_SCOPE)
    endif()
endmacro()
macro (add_executable _name)
    # invoke built-in add_executable
    _add_executable(${ARGV})
    if (TARGET ${_name})
        target_link_libraries(${_name} ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
    endif()
endmacro()

PROJECT(LibMIA)
ENABLE_TESTING()
set (LibMIA_VERSION_MAJOR 0)
set (LibMIA_VERSION_MINOR 5)
set(LibMIA_VERSION ${LibMIA_VERSION_MAJOR}.${LibMIA_VERSION_MINOR})

#set C++11 flag
# GCC specifics
if(CMAKE_COMPILER_IS_GNUCXX OR MINGW OR CYGWIN)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline-functions -ftemplate-depth-128")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -std=gnu++0x -Drestrict_libmia=__restrict__ ")  
  
endif()

# Visual Studio specifics - C++11 is on automatically
if(MSVC)
  SET(INTEL_COMPILER TRUE CACHE BOOL "if you are using Intel compiler, will add intel-specific optimization flags. You still need to manually tell Visual Studio to use Intel compiler if you are using that IDE")
  #msvc doesn't support constexpr yet, so we need to define an empty macro and also add an appropriate compiler flag
  
  if(INTEL_COMPILER)
		
		add_definitions(-Drestrict_libmia=restrict)
		
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Ob2 -Qipo -Zc:forScope -MD -nologo -Qopt-matmul -Ot -Oi -Qansi_alias -QxHost -Quse-intel-optimized-headers")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qrestrict -Qstd=c++11")
		
  else(INTEL_COMPILER)
		add_definitions(-Drestrict_libmia=__restrict)
		add_definitions(-DNO_CONST_EXPR_LIBMIA)
		add_definitions(-DNOMINMAX )
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bigobj")
		message("MSVC compiler used. No support for MIA Functions, but MATLAB wrapper can be built")
  endif()
  add_definitions(-D_ALLOW_KEYWORD_MACROS)
  
  # Warning Level 4 for Debug builds 
  list(APPEND CMAKE_CXX_FLAGS_DEBUG "/W4")
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS)
endif()
#FIND_PACKAGE( OpenMP)
#if(OPENMP_FOUND)
	
#	if(INTEL_COMPILER)
#		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qopenmp")
#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qopenmp")
#	else(INTEL_COMPILER)
#		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#	endif()
		
	
#	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()
SET(USE_MKL FALSE CACHE BOOL "dense tensor multiplication routines can use mkl")
if(USE_MKL)
	SET(MKL_ROOT CACHE PATH /path/to/mkl)
	include_directories("${MKL_ROOT}/include") 
	LINK_DIRECTORIES(${MKL_ROOT}/lib/intel64/)
	set(MKL_LIBS "mkl_intel_lp64")
	set(MKL_LIBS1 "mkl_core")
	if(OPENMP_FOUND)
		set(MKL_LIBS2 "mkl_intel_thread")	
	else(OPENMO_FOUND)
		set(MKL_LIBS2 "mkl_sequential")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_USE_MKL_ALL")
endif()


# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH  "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/LibMIA)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH  "Installation directory for CMake files")

# # Make relative paths absolute (needed later on)
# foreach(p LIB BIN INCLUDE CMAKE)
  # set(var INSTALL_${p}_DIR)
  # if(NOT IS_ABSOLUTE "${${var}}")
    # set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  # endif()
# endforeach()

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmia_doxygen.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
MESSAGE("Doxygen found, use make doc to build class api documentation")

endif(DOXYGEN_FOUND)



mark_as_advanced(CLEAR BOOST_ROOT)

add_subdirectory (src)


add_library(LibMIAHeaders INTERFACE)
target_include_directories(LibMIAHeaders INTERFACE
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR}/src/multi_index_array>
  $<INSTALL_INTERFACE:include/LibMIAHeaders>
)

install(TARGETS LibMIAHeaders EXPORT libMIAExport)
install(EXPORT libMIAExport NAMESPACE Upstream::
  DESTINATION lib/cmake/LibMIA
)

set(LibMIA_headers ${PROJECT_BINARY_DIR}/MIAConfig.h ${SRCS})
add_custom_target(LibMIA SOURCES ${LibMIA_headers})

install(FILES
    ${LibMIA_headers}
  DESTINATION include/LibMIAHeaders
)

# # Add all targets to the build-tree export set
# export(TARGETS LibMIA
  # FILE "${PROJECT_BINARY_DIR}/LibMIALibraryDepends.cmake")
  
# # Export the package for use from the build-tree
# # (this registers the build-tree with a global CMake-registry)
# export(PACKAGE LibMIA)

# # Create a LibMIABuildTreeSettings.cmake file for the use from the build tree
# configure_file(LibMIABuildTreeSettings.cmake.in
  # "${PROJECT_BINARY_DIR}/LibMIABuildTreeSettings.cmake" @ONLY)
 
# # Create the LibMIAConfig.cmake and LibMIAConfigVersion files
# file(RELATIVE_PATH CONF_REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   # "${INSTALL_INCLUDE_DIR}")
# configure_file(LibMIAConfig.cmake.in
  # "${PROJECT_BINARY_DIR}/LibMIAConfig.cmake" @ONLY)
# configure_file(LibMIAConfigVersion.cmake.in
  # "${PROJECT_BINARY_DIR}/LibMIAConfigVersion.cmake" @ONLY)
 
# # Install the LibMIAConfig.cmake and LibMIAConfigVersion.cmake
# install(FILES
  # "${PROJECT_BINARY_DIR}/LibMIAConfig.cmake"
  # "${PROJECT_BINARY_DIR}/LibMIAConfigVersion.cmake"
  # DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# # Install the export set for use with the install-tree
# install(EXPORT LibMIALibraryDepends DESTINATION
  # "${INSTALL_CMAKE_DIR}" )
