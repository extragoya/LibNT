function C=plus(A,B)

A.inds=A.inds-1;
B.inds=B.inds-1;
[c_data c_inds]=sparselattice_merge(A.vals,A.inds,[A.m A.n A.p],B.vals,B.inds,[B.m B.n B.p],0);
C=SparseLattice(c_inds,c_data,A.m,A.n,A.p);
A.inds=A.inds+1;
B.inds=B.inds+1;




