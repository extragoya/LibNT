function [C, c_data]=dense_mia_merge_test_work(a_dims,a_indices,b_indices,permute_idx,doplus)
%a_indices, b_indices, the string of indices, e.g. 'ijk'. permute_idx, the
%actual numerical indices by which A and B are matched. That is
%A{1,2,3}=B{permute_idx(1),permute_idx(2),permute_idx(3)};


a_data=rand(a_dims);
inverse_idx=permute_idx;
inverse_idx( permute_idx ) = 1:length(permute_idx);
b_dims=a_dims(inverse_idx);
b_data=rand(b_dims);
A=MIA(a_data);
B=MIA(b_data);
%testing assignment with indices has its own dedicated test, so just pull
%the MIA from resulting Expr class
if(doplus)
    CExpr=A(a_indices)+B(b_indices);
else
   CExpr=A(a_indices)-B(b_indices); 
end
C=CExpr.m_mia;
b_data=permute(b_data,permute_idx);
if(doplus)
    c_data=a_data+b_data;
else
    c_data=a_data-b_data;
end




