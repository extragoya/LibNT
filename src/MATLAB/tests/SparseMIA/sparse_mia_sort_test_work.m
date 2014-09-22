function [A, A2]=sparse_mia_sort_test_work(dims)


a_data=rand(dims);
a_data(a_data<0.5)=0;

DenseA=MIA(a_data);
A=SparseMIA(DenseA);
A2=SparseMIA(DenseA);
P = randperm(A.nnz);
A.indices=A.indices(P);
A2.indices=A2.indices(P);
A.isSorted=false;
A2.isSorted=false;
A.sort();
%verify against MATLAB's sort
if(A2.nnz)
    [A2.indices, idx]=sort(A2.indices);
    A2.data=A2.data(idx);
end

   











