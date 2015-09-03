function [C c_data]=dense_lattice_solve_test_work(m1,n1,m2,n2,p)

a_data=rand(m1,n1,p);
b_data=rand(m2,n2,p);

A=Lattice(a_data);
B=Lattice(b_data);

C=A\B;

c_data=zeros(n1,n2,p);

for i=1:p
    c_data(:,:,i)=a_data(:,:,i)\b_data(:,:,i);
end





