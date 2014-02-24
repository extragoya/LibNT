function Expr=make_expr(obj, indices)

split_indices=MIA.make_cell_indices(indices);
if length(split_indices)~=obj.order
       error('Number of indices must be the same as MIA order') 
end
Expr=MIAExpr(obj,split_indices);
