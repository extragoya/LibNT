function C=mtimes(A,B)

if A.n~=B.m
    error('Number of columns in A must equal number of rows in B')
elseif A.p~=B.p
    error('Depth of both lattices must be identical')
end
c_data=DenseLatticeMultMex(A.vals,B.vals);
C=Lattice(c_data);