function [C, DenseC]=sparse_mia_mult_test_work(a_dims,b_dims,a_indices,b_indices,a_lexOrder,b_lexOrder,a_sorted,b_sorted)

a_data=rand(a_dims);
a_data(a_data<0.5)=0;
DenseA=MIA(a_data);
A=SparseMIA(DenseA);
A=A.permute(a_lexOrder);
if(~a_sorted)
    P = randperm(A.nnz);
    A.indices=A.indices(P);
    A.data=A.data(P);
    A.isSorted=false;
end
b_data=rand(b_dims);
b_data(b_data<0.5)=0;
DenseB=MIA(b_data);
B=SparseMIA(DenseB);
B=B.permute(b_lexOrder);
if(~b_sorted)
    P = randperm(B.nnz);
    B.indices=B.indices(P);
    B.data=B.data(P);
    B.isSorted=false;
end

%we will test indexed assignment in its own test, so just save the
%resulting MIAExpr and pull its MIA.
CExpr=A(a_indices)*B(b_indices);
C=CExpr.m_mia;
CExpr=DenseA(a_indices)*DenseB(b_indices);
DenseC=CExpr.m_mia;




