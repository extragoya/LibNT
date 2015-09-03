function C= do_minus(A,B,permute_idx)
%perform addition based on the permutation idx of B
if ~isa(A,'DenseNT') || ~isa(B,'DenseNT')
    error('SparseNTs can only be substracted with NT classes');
end

if isa(A,'SparseNT') && isa(B,'SparseNT')
    
    C=do_merge(A,B,permute_idx,1);
else
    error('Sparse NT plus dense NT functionality not supported yet');
end

end

