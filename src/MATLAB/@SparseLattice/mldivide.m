function B=mldivide(A,C)

if A.p~=C.p
    error('Depth of both lattices must be identical')
end

if isa(C,'SparseLattice')
    A.inds=A.inds-1;
    C.inds=C.inds-1;
    b_vals=SparseLatticeSolveMex(A.vals,A.inds,[A.m A.n A.p],C.vals,C.inds,[C.m C.n C.p]);
    
    B=Lattice(b_vals);
    A.inds=A.inds+1;
    C.inds=C.inds+1;
end
