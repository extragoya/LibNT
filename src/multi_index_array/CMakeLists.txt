cmake_minimum_required (VERSION 2.6)
set(LibMIA_headers
     DenseLattice.h
     DenseLatticeBase.h
     DenseMIA.h
     DenseMIABase.h
	 ExprUtil.h
     Lattice.h
     Index.h
	 IndexUtil.h
     LibMIAAlgorithm.h
	 LibMIAException.h
     MappedDenseLattice.h
     MappedSparseLattice.h
     MIA.h
     MIA_Expr.h	 
     Packer.h
	 PermuteIterator.h
     SparseLattice.h
     SparseLatticeBase.h
	 SparseMIA.h
	 SparseMIABase.h
     tupleit.hh
     UniqueLabel.h
     Util.h  
	${PROJECT_BINARY_DIR}/MIAConfig.h	 
     )
add_library(LibMIA STATIC ${LibMIA_headers})
set_target_properties(LibMIA PROPERTIES LINKER_LANGUAGE CXX)
if(LIBMIA_USE_SPARSE_SOLVE)
	TARGET_LINK_LIBRARIES(LibMIA ${SuperLU_LIB} ${blas_LIB})
endif(LIBMIA_USE_SPARSE_SOLVE)

set_target_properties(LibMIA PROPERTIES
  PUBLIC_HEADER "${LibMIA_headers};${CMAKE_BINARY_DIR}/MIAConfig.h")
 
install(TARGETS LibMIA
  # IMPORTANT: Add the LibMIA library to the "export-set"
  EXPORT LibMIALibraryDepends  
  #we make a dummy LibMIA library, so another library can link to it and have it resolve its own library dependencies.
  ARCHIVE DESTINATION "${INSTALL_INCLUDE_DIR}" OPTIONAL
  #LIBRARY DESTINATION "${INSTALL_LIB_DIR}/shlib" 
  #ARCHIVE DESTINATION include OPTIONAL 
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
    )

