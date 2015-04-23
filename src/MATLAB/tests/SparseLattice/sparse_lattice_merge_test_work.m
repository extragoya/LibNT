function [C, DenseC]=sparse_lattice_merge_test_work(m,n,p,op)


a_data=rand(m,n,p);
b_data=rand(m,n,p);
a_data(a_data<0.5)=0;
b_data(b_data<0.5)=0;
DenseA=Lattice(a_data);
DenseB=Lattice(b_data);
A=SparseLattice(a_data);
B=SparseLattice(b_data);
if(op==1)
    DenseC=DenseA+DenseB;
    C=A+B;    
else
    DenseC=DenseA-DenseB;
    C=A-B;
end
   











