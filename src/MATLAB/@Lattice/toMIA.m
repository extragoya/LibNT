function A_MIA=toMIA(A,row_size,col_size,depth_size)
%converts a lattice to MIA based by specifying new sizes of row, columns, 
% and depth ranges

% 
new_dims=[row_size col_size depth_size];
if(find(new_dims==1))
    error('Cannot input singleton dimensions');
end

A_data=squeeze(A.vals);
if length(new_dims)>1
    A_data=reshape(A_data,new_dims);
end
A_MIA=DenseNT(A_data,A.solveInfo);