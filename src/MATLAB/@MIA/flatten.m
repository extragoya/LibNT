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

row_dims=dims(row_idx);
col_dims=dims(col_idx);
row_length=prod(row_dims);
col_length=prod(col_dims);

error_check_flatten(A,row_idx,col_idx);



B=A.toLattice(row_idx, col_idx, [],[row_length col_length 1]);
B=squeeze(B.vals);
