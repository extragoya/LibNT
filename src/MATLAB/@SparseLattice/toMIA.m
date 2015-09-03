function A_MIA=toMIA(A,row_size,col_size,depth_size,transposed)
%converts a lattice to MIA based by specifying new sizes of row, columns, 
% and depth ranges. Does not perform a permute. However, if transposed is
% set to true, it will assume lat is in RowMajor form

% 
if nargin==4
    transposed=false;
end
if(~transposed)
    new_dims=[row_size col_size depth_size];


    A_MIA=SparseNT(A.vals,A.inds,new_dims,1:length(new_dims),true);
else
    new_dims=[col_size row_size  depth_size];
    
    startLinIdx=length(col_size)+1;
    startDepthIdx=startLinIdx+length(row_size);
    A_MIA=SparseNT(A.vals,A.inds,new_dims,[startLinIdx:startLinIdx+length(row_size)-1 1:startLinIdx-1 startDepthIdx:length(new_dims)],true);
end

