function B=toLattice(A,row_idx,col_idx,depth_idx)
a_size=size(A);
a_row_size=prod(a_size(row_idx));
a_col_size=prod(a_size(col_idx));
a_depth_size=prod(a_size(depth_idx));

A_data=A.data;
%extract data values and reorder elements 
if numel([row_idx col_idx depth_idx])>1
    A_data=permute(A_data,[row_idx col_idx depth_idx]);
end

%flatten reordered elements so that their indices correspond to appropriate
%rows, columns, or tabs of lattices
A_data=reshape(A_data,[a_row_size a_col_size a_depth_size]);
B=Lattice(A_data);