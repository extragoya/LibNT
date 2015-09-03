cmake_minimum_required (VERSION 2.6)
PROJECT(LibMIA_example)

#first do some compiler checkes as each compiler has different flags. So far, this cmake file should set you up for GCC, MSVC, and Intel. Other compilers can be added
#turn on C++11 flag
if(CMAKE_COMPILER_IS_GNUCXX OR MINGW OR CYGWIN)
	#gcc accepts a tag for restricted pointers (https://gcc.gnu.org/onlinedocs/gcc/Restricted-Pointers.html), other compilers have their own flags and rules, so you can adjust accordingly.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x ${CMAKE_CXX_FLAGS} -Drestrict_libmia=__restrict__")  
  
  #You can also add your own more aggressive optimizations here if you wish.
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline-functions -ftemplate-depth-128") 
endif()

# Visual Studio specifics - C++11 is on automatically
if(MSVC)
  #Cmake doesn't have a good way to detect when you're using the Intel compiler with microsoft visual studio, so this is a manual check box
  SET(INTEL_COMPILER TRUE CACHE BOOL "if you are using Intel compiler, will add intel-specific optimization flags. You still need to manually tell Visual Studio to use Intel compiler if you are using that IDE")  
  if(INTEL_COMPILER)		
		add_definitions(-Drestrict_libmia=restrict)		
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Ob2 -Qipo -Zc:forScope -MD -nologo -Qopt-matmul -Ot -Oi -Qansi_alias -QxHost -Quse-intel-optimized-headers") #full optimization
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qrestrict -Qstd=c++11")		
  else(INTEL_COMPILER)		
		add_definitions(-Drestrict_libmia=__restrict)		#more aggressive optimizations can be added here		
  endif()
  add_definitions(-D_ALLOW_KEYWORD_MACROS)
  
  # Warning Level 4 for Debug builds 
  list(APPEND CMAKE_CXX_FLAGS_DEBUG "/W4")
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS)
endif()

#Eigen (which LibNT relies on) has openmp support, so may as well turn it on if openmp is supported
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)	
	if(INTEL_COMPILER) #cmake's openmp flags will be configured for msvc, so manually change them
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qopenmp")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qopenmp")
	else(INTEL_COMPILER)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()	
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#finally, Eigen also has optional intel mkl support, so if you have the intel mkl installed, you can link to it to make Eigen run faster
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
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_USE_MKL_ALL") #tell eigen to use intel mkl libraries
endif()


#point CMake to your boost directory
SET(BOOST_ROOT CACHE PATH /path/to/boost)
find_package(Boost REQUIRED)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS}) 
endif()

#point CMake to your Eigen directory
SET(Eigen_ROOT  CACHE PATH /path/to/eigen)
include_directories(${Eigen_ROOT}) 
#point CMake to your <libmia_root>/src/multi_index_array  directory
SET(LibMIA_ROOT  CACHE PATH /path/to/libmia_root/src/multi_index_array)
include_directories(${LibMIA_ROOT}) 


#create a target executable called example_cmake using the main.cpp source file, and link to the mkl libraries if you are using them
add_executable(example_cmake main.cpp )
TARGET_LINK_LIBRARIES(example_cmake  ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})