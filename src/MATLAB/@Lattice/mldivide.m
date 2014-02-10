function B=mldivide(A,C)

if A.p~=C.p
    error('Depth of both lattices must be identical')
end


[b_vals, res]=DenseLatticeSolveMex(A.vals,C.vals);


B=Lattice(b_vals);
B.solveInfo=res;