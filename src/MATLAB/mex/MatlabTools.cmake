# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  MatlabTools.cmake
# @brief Enables use of MATLAB Compiler and build of MEX-files.
#
# @ingroup CMakeTools
##############################################################################



# ----------------------------------------------------------------------------
## @brief Determine extension of MEX-files for this architecture.
#
# @param [out] ARGN The first argument ARGV0 is set to the extension of
#                   MEX-files (excluding '.'). If the CMake variable MEX_EXT
#                   is set, its value is returned. Otherwise, this function
#                   tries to determine it from the system information.
#                   If the extension could not be determined, an empty string
#                   is returned. If no argument is given, the extension is
#                   cached as the variable MEX_EXT.
#
# @returns Sets the variable named by the first argument to the
#          platform-specific extension of MEX-files.
#
# @ingroup CMakeUtilities
function (basis_mexext)
  # default return value
  set (MEXEXT "${MEX_EXT}")
  # use MEXEXT if possible
  if (NOT MEXEXT AND MATLAB_MEXEXT_EXECUTABLE)
    execute_process (
      COMMAND         "${MATLAB_MEXEXT_EXECUTABLE}"
      RESULT_VARIABLE RETVAL
      OUTPUT_VARIABLE MEXEXT
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (RETVAL)
      set (MEXEXT "")
    endif ()
  endif ()
  # otherwise, determine extension given CMake variables describing the system
  if (NOT MEXEXT)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR
          CMAKE_SIZEOF_VOID_P MATCHES 8)
        set (MEXEXT "mexa64")
      elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86" OR
              CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
        set (MEXEXT "mexglx")
      endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES "Windows")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR
          CMAKE_SIZEOF_VOID_P    MATCHES 8)
        set (MEXEXT "mexw64")
      elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86" OR
              CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
        set (MEXEXT "mexw32")
      endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR
          CMAKE_SIZEOF_VOID_P    MATCHES 8)
        set (MEXEXT "mexmaci64")
      else ()
        set (MEXEXT "mexmaci")
      endif ()
    elseif (CMAKE_SYSTEM_NAME MATCHES "SunOS")
      set (MEXEXT "mexs64")
    endif ()
  endif ()
  # return value
  if (ARGC GREATER 0)
    set ("${ARGV0}" "${MEXEXT}" PARENT_SCOPE)
  else ()
    if (NOT DEFINED MEX_EXT)
      set (MARKIT 1)
    else ()
      set (MARKIT 0)
    endif ()
    set (MEX_EXT "${MEXEXT}" CACHE STRING "The extension of MEX-files for this architecture." FORCE)
    if (MARKIT)
      mark_as_advanced (MEX_EXT)
    endif ()
  endif ()
endfunction ()

