cmake_minimum_required (VERSION 2.6)

add_subdirectory(Lattice)
add_subdirectory(MIA)
if(USE_SPARSE_SOLVE)
  add_subdirectory(SparseLattice)
endif(USE_SPARSE_SOLVE)