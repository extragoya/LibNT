function A = sort(A)
%Sorts the SparseNT. 
%   Sorts the SparseMIA based on its index values. Unlike permute, assumes 
%   indices are unsorted to begin with

if(~A.isSorted && A.nnz>0)
    SparseSortMex(A.data,A.indices);
end


end

