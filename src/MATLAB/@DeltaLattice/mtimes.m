function C= mtimes(A,B)



if(isa(A,'SparseLattice') || isa(B,'SparseLattice'))
    C=sparse_times(A,B);
elseif (isa(A,'Lattice') || isa(B,'Lattice'))
    C=dense_times(A,B);
else
    error('Can only multiply a delta lattice with a dense or sparse lattice')
end
  
