function B=flatten(A,row_idx,col_idx)


dims=size(A);
N=numel(dims);
if N==2
    if dims(2)==1
        N=1;
        dims=dims(1);
    end
end

if nargin==1 ||(isempty(row_idx) && isempty(col_idx))
    row_idx=1:N;
    col_idx=[];
elseif nargin==2
    col_idx=[];
end
error_check_flatten(A,row_idx,col_idx);

row_dims=dims(row_idx);
row_length=prod(row_dims);
if ~isempty(row_dims)
    row_idx=A.extract_indices(row_idx);
     
else
    row_idx=ones(length(A.data),1); 
    
end

col_dims=dims(col_idx);
col_length=prod(col_dims);
if ~isempty(col_dims)
    col_idx=A.extract_indices(col_idx);
    
else
    col_idx=ones(length(A.data),1); 
    
end



B=sparse(double(row_idx),double(col_idx),A.data,row_length,col_length,length(A.data));