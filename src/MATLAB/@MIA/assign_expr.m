function obj=assign_expr(obj,indices,Expr)
%When doing subscripted assignment, i.e., A('ijk')=B('kji'), this function
%is called. Right side has to be an MIAExpr class. Function will check that
%indices match properly, and if so, then the actual assignment will be
%performed
otherMIA=Expr.m_mia;

split_indices=MIA.make_cell_indices(indices);
other_indices=Expr.m_indices;
if(length(split_indices)~=length(other_indices))
    error('Orders of two MIAs must be the same to perform assignment');
end
cur_order=length(split_indices);
reg='[a-z]';
assign_order=zeros(cur_order,1);
%try to match every index on LHS
for i=1:cur_order
    cur_a_index=split_indices{i};
    %remove any exclamation points    
    cur_a_index=regexp(cur_a_index,reg,'match'); %unless make_cell_indices screwed up, there should only be one char
    cur_a_index=cur_a_index{1};
    for j=1:cur_order
        cur_b_index=other_indices{j};
        cur_b_index=regexp(cur_b_index,reg,'match');
        cur_b_index=cur_b_index{1};
        if(strcmp(cur_a_index,cur_b_index))
            assign_order(i)=j;
        end
    end
    %if we didn't find a match, throw an error
    if(assign_order(i)==0)
        error(['Index ' cur_a_index ' in LHS was not matched'])
    end

end
obj=assign(obj,otherMIA,assign_order); %perform actual assignment