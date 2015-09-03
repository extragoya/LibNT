function [C, c_data]=dense_mia_mult_test_work(a_dims,b_dims,a_indices,b_indices,a_outer_idx,a_inner_idx,a_inter_idx,b_outer_idx,b_inner_idx,b_inter_idx)
%a_indices, b_indices are hte string indices, e.g. 'ijk'. a_outer_idx,
%a_inter_idx, etc are the actual numerical index values that are undergoing
%each of the three possible mulitplications
a_data=rand(a_dims);
b_data=rand(b_dims);

A=DenseNT(a_data);
B=DenseNT(b_data);

%we will test indexed assignment in its own test, so just save the
%resulting MIAExpr and pull its MIA.
CExpr=A(a_indices)*B(b_indices);
C=CExpr.m_mia;
a_data=permute(a_data,[a_outer_idx a_inner_idx a_inter_idx]);
b_data=permute(b_data,[b_inner_idx b_outer_idx b_inter_idx]);
a_row_size=prod(a_dims(a_outer_idx));
a_col_size=prod(a_dims(a_inner_idx));
a_depth_size=prod(a_dims(a_inter_idx));
b_row_size=prod(b_dims(b_inner_idx));
b_col_size=prod(b_dims(b_outer_idx));
b_depth_size=prod(b_dims(b_inter_idx));
a_data=reshape(a_data,[a_row_size a_col_size a_depth_size]);
b_data=reshape(b_data,[b_row_size b_col_size b_depth_size]);
c_data=zeros(a_row_size,b_col_size,a_depth_size);

for i=1:a_depth_size
    c_data(:,:,i)=a_data(:,:,i)*b_data(:,:,i);
end
new_dims=[a_dims(a_outer_idx) b_dims(b_outer_idx) a_dims(a_inter_idx)];
if (length(new_dims)==1)
    new_dims=[new_dims 1];
end
c_data=reshape(c_data,new_dims);




