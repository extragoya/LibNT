function A_MIA=toMIA(A,row_size,col_size,depth_size)
%converts a lattice to MIA based by specifying new sizes of row, columns, 
% and depth ranges

% 
A_data=A.vals;

A_data=squeeze(A_data);
new_dims=[row_size col_size depth_size];
if length(new_dims)>1
    A_data=reshape(A_data,new_dims);
end
A_data=squeeze(A_data);
A_MIA=MIA(A_data,A.solveInfo);