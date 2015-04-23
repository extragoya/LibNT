function [B, B_test]=dense_mia_solve_test_work(a_dims,b_dims,a_indices,b_indices,LSQR)
%a_indices, b_indices are hte string indices, e.g. 'ijk'. a_outer_idx,
%a_inter_idx, etc are the actual numerical index values that are undergoing
%each of the three possible mulitplications
a_data=rand(a_dims);
b_data=rand(b_dims);

A=MIA(a_data);
B=MIA(b_data);

CExpr=A(a_indices)\B(b_indices); %solve
C=CExpr.m_mia;
%now redo the A*C multiplication to verify C actually solved the equations

if(~LSQR)
    B_test=MIA;
    B_test(b_indices)=A(a_indices)*C(CExpr.m_indices);
else
    D=MIA;
    D(CExpr.m_indices)=A(a_indices)*B(b_indices);
    D2=MIA;
    D2(CExpr.m_indices)=A(a_indices)*(A(a_indices)*C(CExpr.m_indices));
    B=D;
    B_test=D2;
end






