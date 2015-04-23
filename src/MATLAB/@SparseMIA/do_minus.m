function C= do_minus(A,B,permute_idx)
%perform addition based on the permutation idx of B
if ~isa(A,'MIA') || ~isa(B,'MIA')
    error('SparseMIAs can only be multipled with classes or subclasses of MIA');
end

if isa(A,'SparseMIA') && isa(B,'SparseMIA')
    
    C=do_merge(A,B,permute_idx,1);
end

end

