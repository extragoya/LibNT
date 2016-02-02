function B=mldivide(A,C)
B=Lattice;
if A.p~=C.p
    
    error('Depth of both lattices must be identical')
end

if isa(C,'SparseLattice')
    A.inds=A.inds-1;
    C.inds=C.inds-1;
    b_data=SparseLatticeSolveMex(A.vals,A.inds,[A.m A.n A.p],C.vals,C.inds,[C.m C.n C.p]);
    b_data=reshape(b_data,[A.n C.n A.p]);
    B=Lattice(b_data);    
    A.inds=A.inds+1;
    C.inds=C.inds+1;
elseif isa(C,'Lattice') %to do, alter mex function so it can accept dense data
      
    A.inds=A.inds-1;    
    b_data=SparseDenseLatticeSolveMex(A.vals,A.inds,[A.m A.n A.p],C.vals,0);
    b_data=reshape(b_data,[A.n C.n A.p]);
    B=Lattice(b_data);    
    A.inds=A.inds+1;
    
end
