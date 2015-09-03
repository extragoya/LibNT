function permute_idx=pull_merge_indices(A,B)
%Given two NTExpr's will pull how to permute B to get its indices to match
%up with A. E.g., if a_indices is 'ijk' and b_indices is 'jki', will
%produce array [3 1 2].
if ~isa(B,'NTExpr') 
    error('NTExpr''s can only be merged with other NTExpr''s');
end
a_order=A.m_mia.order;
b_order=B.m_mia.order;
if(a_order~=b_order)
    error('Only NTs of the same order can be added');
end

a_indices=A.m_indices;
b_indices=B.m_indices;


permute_idx=zeros(a_order,1);
%try to match every index on LHS
for i=1:a_order
    cur_a_index=a_indices{i};
    %remove any exclamation points    
    letter_start=regexp(cur_a_index,'[a-z]', 'once');
    cur_a_index=cur_a_index(letter_start:end); %unless make_cell_indices screwed up, there should only be one char   
    for j=1:a_order
        cur_b_index=b_indices{j};
        letter_start=regexp(cur_b_index,'[a-z]', 'once');
        cur_b_index=cur_b_index(letter_start:end);        
        if(strcmp(cur_a_index,cur_b_index))
            permute_idx(i)=j;
        end
    end
    %if we didn't find a match, throw an error
    if(permute_idx(i)==0)
        error(['Index ' cur_a_index ' in LHS was not matched'])
    end

end