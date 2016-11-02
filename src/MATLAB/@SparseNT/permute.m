function B = permute(A,newLinIdx)
%Permute the SparseNT. 
%   Permutes the SparseNT so that its linear indices are calculated and
%   sorted in the new lexicographical order. If A is not already sorted, it
%   will just call sort


B=A;
B.data=B.data+1;
B.data=B.data-1;
if(~B.isSorted)
    B=B.changeLexOrder(newLinIdx);
    SparseSortMex(B.data,B.indices);
elseif(~isequal(newLinIdx,B.lexOrder))
    oldLexOrder=B.lexOrder;
    B=B.changeLexOrder(newLinIdx);
    B.indices=B.indices-1;
    SparsePermuteMex(B.data,B.indices,uint8(oldLexOrder),uint8(newLinIdx),int64(B.dims));    
    B.indices=int64(B.indices+1);
    
    B.isSorted=true;
end


end

