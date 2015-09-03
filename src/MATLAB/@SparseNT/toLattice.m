function B=toLattice(A,row_idx,col_idx,depth_idx)

a_size=size(A);
a_row_size=prod(a_size(row_idx));
a_col_size=prod(a_size(col_idx));
a_depth_size=prod(a_size(depth_idx));


%extract data values and reorder elements 
if numel([row_idx col_idx depth_idx])>1
    A=A.permute([row_idx col_idx depth_idx]);
end

B=SparseLattice(A.indices,A.data,a_row_size,a_col_size,a_depth_size,true);