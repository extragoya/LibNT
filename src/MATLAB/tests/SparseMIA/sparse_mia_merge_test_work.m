function [C, DenseC]=sparse_mia_merge_test_work(a_dims,a_indices,a_lexOrder,a_sorted,b_dims,b_indices,b_lexOrder,b_sorted,do_plus)
%a_indices, b_indices, the string of indices, e.g. 'ijk'. permute_idx, the
%actual numerical indices by which A and B are matched. That is
%A{1,2,3}=B{permute_idx(1),permute_idx(2),permute_idx(3)};

a_data=rand(a_dims);
a_data(a_data<0.5)=0;
DenseA=DenseNT(a_data);
A=SparseNT(DenseA);
A=A.permute(a_lexOrder);
if(~a_sorted)
    P = randperm(A.nnz);
    A.indices=A.indices(P);
    A.data=A.data(P);
    A.isSorted=false;
end
b_data=rand(b_dims);
b_data(b_data<0.5)=0;
DenseB=DenseNT(b_data);
B=SparseNT(DenseB);
B=B.permute(b_lexOrder);
if(~b_sorted)
    P = randperm(B.nnz);
    B.indices=B.indices(P);
    B.data=B.data(P);
    B.isSorted=false;
end
%testing assignment with indices has its own dedicated test, so just pull
%the MIA from resulting Expr class
if(do_plus)
    CExpr=A(a_indices)+B(b_indices);
    CExpr2=DenseA(a_indices)+DenseB(b_indices);
else
    CExpr=A(a_indices)-B(b_indices); 
    CExpr2=DenseA(a_indices)-DenseB(b_indices);

end
C=CExpr.m_mia;
DenseC=CExpr2.m_mia;




