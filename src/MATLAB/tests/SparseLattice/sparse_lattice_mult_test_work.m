function [C DenseC]=sparse_lattice_mult_test_work(m1,n1,m2,n2,p)

a_data=rand(m1,n1,p);
b_data=rand(m2,n2,p);
a_data(a_data<0.5)=0;
b_data(b_data<0.5)=0;
DenseA=Lattice(a_data);
DenseB=Lattice(b_data);
A=SparseLattice(a_data);
B=SparseLattice(b_data);

C=A*B;
DenseC=DenseA*DenseB;






