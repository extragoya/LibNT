cmake_minimum_required (VERSION 2.6)

ENABLE_TESTING()
add_subdirectory(Lattice)
add_subdirectory(MIA)
add_subdirectory(DenseMIA)
#this should be changed 
if(USE_SPARSE_SOLVE)
  add_subdirectory(SparseLattice)
endif(USE_SPARSE_SOLVE)
