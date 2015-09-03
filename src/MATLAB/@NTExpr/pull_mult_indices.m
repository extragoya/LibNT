function [a_inner_idx, a_inter_idx, a_outer_idx, b_inner_idx, b_inter_idx, b_outer_idx]=pull_mult_indices(A,B)
%Given two NTExpr's will pull the sequence of inner, inter and outer
%product indices expressed by the two expressions
if ~isa(B,'NTExpr') 
    error('NTExpr''s can only be multipled with other NTExpr''s');
end



%find common indices
%unfortunately it's faster to use set intersect's default option, and then
%manually sort arrays later based on the order of indices in a
a_all=1:length(A.m_indices);
b_all=1:length(B.m_indices);
b_no_intersect= true(1,length(B.m_indices));
[a_idx, b_idx]=ismember(A.m_indices,B.m_indices);
common=A.m_indices(a_idx);
common_a_idx=a_all(a_idx);
common_b_idx=b_idx(a_idx);
b_no_intersect(common_b_idx)=0;
%find common indices with !
found_entrywise = regexp(common,'!');
found_entrywise=cellfun('isempty',found_entrywise);

%indices that have exclamation point are entry-wise products
a_inter_idx=common_a_idx(~found_entrywise);
b_inter_idx=common_b_idx(~found_entrywise);
%indices that don't are inner products
a_inner_idx=common_a_idx(found_entrywise);
b_inner_idx=common_b_idx(found_entrywise);
%get outerproduct indices of both sets
a_outer_idx=a_all(~a_idx);
outer_a=A.m_indices(a_outer_idx);
b_outer_idx=b_all(b_no_intersect);
outer_b=B.m_indices(b_outer_idx);
%now remove exclamation points, to do error checking

only_letter2= regexprep(outer_a,'!','');
only_letter3= regexprep(outer_b,'!','');
%if any index was matched between the tow outer product sets, stripped of exclamation marks, then the user
%must have mismathed exclamation points, e.g., i and !i
test=ismember(only_letter2,only_letter3);
if(sum(test)>0)
  error('Indices must match with both their letter characters and number of exclamation points.')
end

end