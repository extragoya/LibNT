function Expr=make_expr(obj, indices)
if (ischar(indices))
    split_indices=DenseNT.make_cell_indices(indices);
    
elseif iscell(indices)
    split_indices=indices;
end
if length(split_indices)~=obj.order
    error('Number of indices must be the same as MIA order') 
end
Expr=NTExpr(obj,split_indices);