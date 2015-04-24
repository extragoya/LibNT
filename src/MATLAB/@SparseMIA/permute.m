function A = permute(A,newLinIdx)
%Permute the SparseMIA. 
%   Permutes the SparseMIA so that its linear indices are calculated and
%   sorted in the new lexicographical order. If A is not already sorted, it
%   will just call sort



if(~A.isSorted)
    A=A.changeLexOrder(newLinIdx);
    SparseSortMex(A.data,A.indices);
elseif(~isequal(newLinIdx,A.lexOrder))
    oldLexOrder=A.lexOrder;
    A=A.changeLexOrder(newLinIdx);
    A.indices=A.indices-1;
    SparsePermuteMex(A.data,A.indices,uint8(oldLexOrder),uint8(newLinIdx),int64(A.dims));    
    A.indices=int64(A.indices+1);
    A.isSorted=true;
end


end

