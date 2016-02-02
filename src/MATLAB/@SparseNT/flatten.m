function matrix = flatten(obj,row_inds,col_inds)
%FLATTEN flattens MIA into matrix. 
%   When supplied with indices that correponds to rows and columns, will
%   turn MIA data into 2D form. If only row_inds or col_inds are supplied,
%   then will return vector form

dims=size(obj);
N=numel(dims);
if N==2
    if dims(2)==1
        N=1;        
    end
end

if nargin==1 ||(isempty(row_inds) && isempty(col_inds))
    row_inds=1:N;
    col_inds=[];
elseif nargin==2
    col_inds=[];
end


all_indices=[row_inds col_inds];
all_indices_sort=sort(all_indices);
if(~isequal(all_indices_sort,[1:obj.order]))
    error('You must designate all indices, from 1:degree, as being in either row or column indices');
end

matrix=obj.toLattice(row_inds, col_inds, []);
j=idivide(matrix.inds-1,matrix.m)+1;
i=matrix.inds-(j-1)*matrix.m;
matrix = sparse(double(i),double(j),matrix.vals,double(matrix.m),double(matrix.n)) ;