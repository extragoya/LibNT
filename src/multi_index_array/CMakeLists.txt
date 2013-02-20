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
     LibMIAException.h
     MappedDenseLattice.h
     MappedSparseLattice.h
     MIA.h
     MIA_Expr.h
     Packer.h
     SparseLattice.h
     SparseLatticeBase.h
     tupleit.hh
     UniqueLabel.h
     Util.h
     
     )
add_library(LibMIA STATIC ${LibMIA_headers})
set_target_properties(LibMIA PROPERTIES LINKER_LANGUAGE CXX)
if(USE_SPARSE_SOLVE)
	TARGET_LINK_LIBRARIES(LibMIA ${SuperLU_LIB} ${blas_LIB})
endif(USE_SPARSE_SOLVE)

set_target_properties(LibMIA PROPERTIES
  PUBLIC_HEADER "${LibMIA_headers};${CMAKE_BINARY_DIR}/MIAConfig.h")
 
install(TARGETS LibMIA
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT LibMIALibraryDepends  
  ARCHIVE DESTINATION "${INSTALL_INCLUDE_DIR}" OPTIONAL
  #LIBRARY DESTINATION "${INSTALL_LIB_DIR}/shlib" 
  #ARCHIVE DESTINATION include OPTIONAL 
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
    )
#add_custom_target(LibMIA SOURCES ${LibMIA_headers})
