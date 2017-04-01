cmake_minimum_required (VERSION 2.8)
include (GenerateExportHeader)
include (MatlabTools.cmake)
#get the correct extension
basis_mexext()
#makes sure symbols aren't exported by default (probably not an issue when creating mex files, since there's not many symbols anyway)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)



#MSVC has a warning when you use std::copy with pointers. Can't see a way to avoid this right now with the mex interface, so turn it off here.
if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_SCL_SECURE_NO_WARNINGS")
endif(MSVC)
set (SPARSE_CHECK_SOURCES check_sparse.cpp check_sparse.h)
set (DENSE_LATTICE_MULT_SOURCES denselattice_mul.cpp)
set (DENSE_LATTICE_SOLVE_SOURCES denselattice_solve.cpp)
set (SPARSE_LATTICE_MULT_SOURCES sparselattice_mul.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_LATTICE_SOLVE_SOURCES sparselattice_solve.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_MERGE_SOURCES sparse_merge.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_SORT_SOURCES sparse_sort.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_LINIDX_SOURCES sparse_change_linIdx.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_PERMUTE_SOURCES sparse_permute.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_DENSE_LATTICE_MULT_SOURCES sparse_denselattice_mul.cpp ${SPARSE_CHECK_SOURCES})
set (SPARSE_DENSE_LATTICE_SOLVE_SOURCES sparse_denselattice_solve.cpp ${SPARSE_CHECK_SOURCES})
LINK_DIRECTORIES(${MATLAB_LIB_DIR})
add_library(libmx SHARED IMPORTED)
set_property(TARGET libmx PROPERTY IMPORTED_LOCATION "${MATLAB_LIB_DIR}/libmx.so")
add_library(libmex SHARED IMPORTED)
set_property(TARGET libmex PROPERTY IMPORTED_LOCATION "${MATLAB_LIB_DIR}/libmex.so")
add_library(libmat SHARED IMPORTED)
set_property(TARGET libmat PROPERTY IMPORTED_LOCATION "${MATLAB_LIB_DIR}/libmat.so")
#find_library(libmx_location libmx ${MATLAB_LIB_DIR})
#message("libmx is at: ${libmx_location}")

#we set the EXPORT_MACRO_NAME to DLL_EXPORT_SYM, b/c mex.h uses that define. This allows us to export dll symbols without .def file
add_library(DenseLatticeMultMex SHARED ${DENSE_LATTICE_MULT_SOURCES})
add_library(DenseLatticeSolveMex SHARED ${DENSE_LATTICE_SOLVE_SOURCES})
add_library(SparseLatticeMultMex SHARED ${SPARSE_LATTICE_MULT_SOURCES})
add_library(SparseLatticeSolveMex SHARED ${SPARSE_LATTICE_SOLVE_SOURCES})
add_library(SparseMergeMex SHARED ${SPARSE_MERGE_SOURCES})
add_library(SparseSortMex SHARED ${SPARSE_SORT_SOURCES})
add_library(SparseChangeLinIdxMex SHARED ${SPARSE_LINIDX_SOURCES})
add_library(SparsePermuteMex SHARED ${SPARSE_PERMUTE_SOURCES})
add_library(SparseDenseLatticeMultMex SHARED ${SPARSE_DENSE_LATTICE_MULT_SOURCES})
add_library(SparseDenseLatticeSolveMex SHARED ${SPARSE_DENSE_LATTICE_SOLVE_SOURCES})
#Generate a header file that creates the proper export macro
#I believe this function needs a target name when you call it, so I supply one, even though I want the MexSymbols_Export.h file to be used by all mex files.
GENERATE_EXPORT_HEADER( DenseLatticeMultMex
             BASE_NAME DenseLatticeMultMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME DenseLatticeMultMex_Export.h
             STATIC_DEFINE DenseLatticeMultMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( DenseLatticeSolveMex
             BASE_NAME DenseLatticeSolveMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME DenseLatticeSolveMex_Export.h
             STATIC_DEFINE DenseLatticeSolveMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseLatticeMultMex
             BASE_NAME SparseLatticeMultMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseLatticeMultMex_Export.h
             STATIC_DEFINE SparseLatticeMultMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseLatticeSolveMex
             BASE_NAME SparseLatticeSolveMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseLatticeSolveMex_Export.h
             STATIC_DEFINE SparseLatticeSolveMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseMergeMex
             BASE_NAME SparseMergeMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseMergeMex_Export.h
             STATIC_DEFINE SparseMergeMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseSortMex
             BASE_NAME SparseSortMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseSortMex_Export.h
             STATIC_DEFINE SparseSortMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseChangeLinIdxMex
             BASE_NAME SparseChangeLinIdxMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseChangeLinIdxMex_Export.h
             STATIC_DEFINE SparseChangeLinIdxMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparsePermuteMex
             BASE_NAME SparsePermuteMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparsePermuteMex_Export.h
             STATIC_DEFINE SparsePermuteMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseDenseLatticeMultMex
             BASE_NAME SparseDenseLatticeMultMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseDenseLatticeMultMex_Export.h
             STATIC_DEFINE SparseDenseLatticeMultMex_BUILT_AS_STATIC
)
GENERATE_EXPORT_HEADER( SparseDenseLatticeSolveMex
             BASE_NAME SparseDenseLatticeSolveMexSymbols
             EXPORT_MACRO_NAME DLL_EXPORT_SYM
             EXPORT_FILE_NAME SparseDenseLatticeSolveMex_Export.h
             STATIC_DEFINE SparseDenseLatticeSolveMex_BUILT_AS_STATIC
)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
TARGET_LINK_LIBRARIES(DenseLatticeMultMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(DenseLatticeSolveMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseLatticeMultMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseLatticeSolveMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseMergeMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseSortMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseChangeLinIdxMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparsePermuteMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseDenseLatticeMultMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
TARGET_LINK_LIBRARIES(SparseDenseLatticeSolveMex libmx libmex libmat ${MKL_LIBS} ${MKL_LIBS1} ${MKL_LIBS2})
#Change the extension to the appropriate MATLAB one
set_target_properties(DenseLatticeMultMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(DenseLatticeSolveMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseLatticeMultMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseLatticeSolveMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseMergeMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseSortMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseChangeLinIdxMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparsePermuteMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseDenseLatticeMultMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")
set_target_properties(SparseDenseLatticeSolveMex PROPERTIES SUFFIX  .${MEX_EXT} PREFIX "")

add_custom_target(MATLAB)
add_dependencies(MATLAB DenseLatticeMultMex DenseLatticeSolveMex SparseLatticeMultMex SparseLatticeSolveMex SparseMergeMex SparseSortMex SparseChangeLinIdxMex SparsePermuteMex SparseDenseLatticeMultMex SparseDenseLatticeSolveMex)

