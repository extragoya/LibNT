function [B, b_data]=dense_mia_assign_test_work(a_dims,a_indices,b_indices,permute_idx)
%a_indices, b_indices, the string of indices, e.g. 'ijk'. permute_idx, the
%actual numerical indices by which A and B are matched. That is
%A{1,2,3}=B{permute_idx(1),permute_idx(2),permute_idx(3)};


a_data=rand(a_dims);


A=MIA(a_data);
B=MIA;
B(b_indices)=A(a_indices);

b_data=permute(a_data,permute_idx);




