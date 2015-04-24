function [C DenseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,dense_first)

a_data=rand(m1,n1,p);
b_data=rand(m2,n2,p);
if dense_first
    A=Lattice(a_data);
    b_data(b_data<0.5)=0;
    B=SparseLattice(b_data);
else
    
    B=Lattice(b_data);
    a_data(a_data<0.5)=0;
    A=SparseLattice(a_data);
end

DenseA=Lattice(a_data);
DenseB=Lattice(b_data);



C=A*B;
DenseC=DenseA*DenseB;






