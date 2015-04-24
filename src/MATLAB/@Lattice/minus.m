function C=minus(A,B)

if A.n~=B.n ||  A.m~=B.m ||  A.p~=B.p
    error('Dimensions of both lattices must be the same')
end
c_data=A.vals-B.vals;
C=Lattice(c_data);