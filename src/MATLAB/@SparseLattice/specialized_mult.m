function C = specialized_mult(A,B,algorithm)
%will execute one of the specialized multiplication algorithms
if algorithm<0 || algorithm>5
    error('Multiplication algorithm can only be between 1 and 5');
end
m=A.m;

p=A.p;
A.inds=A.inds-1;
B.inds=B.inds-1;
B_lat_dims=[B.m B.n B.p];
if algorithm==4
    B_lat_dims=[B.n B.m B.p]; %b/c the Lattice MATLAB classes only support column-major, this is a work-around to input a row-major lattice
end
n=B_lat_dims(2);
[c_data, c_inds]=SparseLatticeMultMex(A.vals,int64(A.inds),[A.m A.n A.p],B.vals,int64(B.inds),B_lat_dims,algorithm);
c_inds=c_inds+1;
C=SparseLattice(c_inds,c_data,m,n,p);
A.inds=A.inds+1;
B.inds=B.inds+1;

end

