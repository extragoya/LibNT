function A_MIA=toMIA(A,row_size,col_size,depth_size)
%converts a lattice to MIA based by specifying new sizes of row, columns, 
% and depth ranges

% 

new_dims=[row_size col_size depth_size];


A_MIA=SparseMIA(A.vals,A.inds,new_dims,1:length(new_dims),true);

