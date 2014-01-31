function B=toLattice(A,row_idx,col_idx,depth_idx,dims)

A_data=A.data;
%extract data values and reorder elements 
if numel([row_idx col_idx depth_idx])>1
    A_data=permute(A_data,[row_idx col_idx depth_idx]);
end

%flatten reordered elements so that their indices correspond to appropriate
%rows, columns, or tabs of lattices
A_data=reshape(A_data,dims);
B=Lattice(A_data);