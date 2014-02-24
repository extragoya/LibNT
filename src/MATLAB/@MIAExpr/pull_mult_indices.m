function [a_inner_idx, a_inter_idx, a_outer_idx, b_inner_idx, b_inter_idx, b_outer_idx]=pull_mult_indices(A,B)
%Given two MIAExpr's will pull the sequence of inner, inter and outer
%product indices expressed by the two expressions
if ~isa(B,'MIAExpr') 
    error('MIAExpr''s can only be multipled with other MIAExpr''s');
end

a_indices=A.m_indices;
b_indices=B.m_indices;
a_outer_idx=[];
b_outer_idx=[];
a_inner_idx=[];
b_inner_idx=[];
a_inter_idx=[];
b_inter_idx=[];

for i=1:length(a_indices)
    a_cur_idx=a_indices{i};
    a_flag=0;
    for j=1:length(b_indices)
        b_cur_idx=b_indices{j};
        if(strcmp(a_cur_idx,b_cur_idx))            
            a_flag=1;
            if(a_cur_idx(1)=='!')
                a_inter_idx=[a_inter_idx i];
                b_inter_idx=[b_inter_idx j];
            else
                a_inner_idx=[a_inner_idx i];
                b_inner_idx=[b_inner_idx j];
            end
            break;
        end
    end
    if(~a_flag)        
        a_outer_idx=[a_outer_idx i];
    end
end
b_used_indices=[b_inner_idx b_inter_idx];
for i=1:length(b_indices)
    if(sum(b_used_indices==i)==0)
        b_outer_idx=[b_outer_idx i];
    end
end