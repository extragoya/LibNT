function [C DenseC]=sparse_mia_solve_test_work(a_dims,b_dims,a_indices,b_indices,a_lexOrder,b_lexOrder,a_sorted,b_sorted)


flag=true;
while(flag)
    a_data=rand(a_dims);
    b_data=rand(b_dims);
    a_data(a_data<0.5)=0;
    b_data(b_data<0.5)=0;
    DenseA=MIA(a_data);
    DenseB=MIA(b_data);    
    CExpr=DenseA(a_indices)\DenseB(b_indices);
    DenseC=CExpr.m_mia;
    if(DenseC.solveInfo>0)
        flag=false; 
    end
end

A=SparseMIA(DenseA);
A=A.permute(a_lexOrder);
if(~a_sorted)
    P = randperm(A.nnz);
    A.indices=A.indices(P);
    A.data=A.data(P);
    A.isSorted=false;
end
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
CExpr=A(a_indices)\B(b_indices);
C=CExpr.m_mia;
CExpr=DenseA(a_indices)\DenseB(b_indices);
DenseC=CExpr.m_mia;






