function matrix = flatten(obj,row_inds,col_inds)
%FLATTEN flattens MIA into matrix. 
%   When supplied with indices that correponds to rows and columns, will
%   turn MIA data into 2D form. If only row_inds or col_inds are supplied,
%   then will return vector form

all_indices=[row_inds col_inds];
all_indices_sort=sort(all_indices);
if(~isequal(all_indices_sort,[1:obj.order]))
    error('You must designate all indices, from 1:degree, as being in either row or column indices');
end

matrix=obj.toLattice(row_inds, col_inds, []);
j=idivide(matrix.inds-1,matrix.m)+1;
i=matrix.inds-(j-1)*matrix.m;
matrix = sparse(double(i),double(j),matrix.vals,double(matrix.m),double(matrix.n)) ;