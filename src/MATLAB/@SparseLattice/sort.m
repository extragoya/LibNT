function A = sort(A)
%Sorts the SparseMIA. 
%   Sorts the SparseMIA based on its index values. Unlike permute, assumes 
%   indices are unsorted to begin with

if(A.nnz>0)
    SparseSortMex(A.vals,A.inds);
end


end

