function B=mldivide(A,C)

if A.p~=C.p
    error('Depth of both lattices must be identical')
end


b_vals=DenseLatticeSolveMex(A.vals,C.vals);


B=Lattice(b_vals);