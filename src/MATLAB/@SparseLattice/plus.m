function C=plus(A,B)

A.inds=A.inds-1;
B.inds=B.inds-1;
[c_data, c_inds]=SparseMergeMex(A.vals,A.inds,B.vals,B.inds,0);
c_inds=c_inds+1;
C=SparseLattice(int64(c_inds),c_data,A.m,A.n,A.p);
A.inds=A.inds+1;
B.inds=B.inds+1;




