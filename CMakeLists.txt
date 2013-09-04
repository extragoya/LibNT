cmake_minimum_required (VERSION 2.6)
PROJECT(LibMIA)
ENABLE_TESTING()
set (LibMIA_VERSION_MAJOR 0)
set (LibMIA_VERSION_MINOR 5)
set(LibMIA_VERSION ${LibMIA_VERSION_MAJOR}.${LibMIA_VERSION_MINOR})

#set C++11 flag
# GCC specifics
if(CMAKE_COMPILER_IS_GNUCXX OR MINGW OR CYGWIN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -std=gnu++0x -Drestrict=__restrict__")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -std=gnu++0x -Drestrict=__restrict__")
  
endif()

# Visual Studio specifics - C++11 is on automatically
if(MSVC)
  # Warning Level 4 for Debug builds
  list(APPEND CMAKE_CXX_FLAGS_DEBUG " /W4")
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_DEBUG)
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

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

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


add_subdirectory (src)

# Add all targets to the build-tree export set
export(TARGETS LibMIA
  FILE "${PROJECT_BINARY_DIR}/LibMIALibraryDepends.cmake")
  
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE LibMIA)

# Create a LibMIABuildTreeSettings.cmake file for the use from the build tree
configure_file(LibMIABuildTreeSettings.cmake.in
  "${PROJECT_BINARY_DIR}/LibMIABuildTreeSettings.cmake" @ONLY)
 
# Create the LibMIAConfig.cmake and LibMIAConfigVersion files
file(RELATIVE_PATH CONF_REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
configure_file(LibMIAConfig.cmake.in
  "${PROJECT_BINARY_DIR}/LibMIAConfig.cmake" @ONLY)
configure_file(LibMIAConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/LibMIAConfigVersion.cmake" @ONLY)
 
# Install the LibMIAConfig.cmake and LibMIAConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/LibMIAConfig.cmake"
  "${PROJECT_BINARY_DIR}/LibMIAConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT LibMIALibraryDepends DESTINATION
  "${INSTALL_CMAKE_DIR}" )