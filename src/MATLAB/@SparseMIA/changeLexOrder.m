function A= changeLexOrder(A,newLexOrder)
% Changes the lexicographical order used to compute linear indices
if(~isequal(A.lexOrder,newLexOrder))
    A.indices=A.indices-1;
    SparseChangeLinIdxMex(A.indices,uint8(A.lexOrder),uint8(newLexOrder),int64(A.dims));
    A.lexOrder=newLexOrder;
    A.indices=A.indices+1;
    A.isSorted=false;
end

end

