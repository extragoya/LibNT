cmake_minimum_required (VERSION 2.6)
set(LibMIA_headers
     DenseLattice.h
     DenseLatticeBase.h
     DenseMIA.h
     DenseMIABase
     Lattice.h
     Index.h
     LatticeException.h
     MappedDenseLattice.h
     MappedSparseLattice.h
     MIA.h
     MIA_Expr.h
     SparseLattice.h
     SparseLatticeBase.h
     tupleit.hh
     UniqueLabel.h
     Util.h
     
     )

add_custom_target(LibMIA SOURCES ${LibMIA_headers})