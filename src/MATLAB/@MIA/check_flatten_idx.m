function check_flatten_idx(A,row_idx,col_idx)

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

if ~isempty(find(row_idx>N, 1)) || ~isempty(find(row_idx<1, 1))
    error(['Row indices must be within 1 and ' num2str(N)]);
end
if ~isempty(find(col_idx>N,1)) || ~isempty(find(col_idx<1, 1))
    error(['Column indices must be within 1 and ' num2str(N)]);
end

if intersect(row_idx,col_idx)
    error('Row and column indices must not overlap')
end

t_dims=ones(N,1);
t_dims(row_idx)=0;
t_dims(col_idx)=0;

if find(t_dims)
    error('Not all indices were set')
end