function Expr=make_expr(obj, indices)

split_indices=make_cell_indices(obj,indices);
Expr=MIAExpr(obj,split_indices);
