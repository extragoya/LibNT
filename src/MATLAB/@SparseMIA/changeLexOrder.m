function A= changeLexOrder(A,newLexOrder)
% Changes the lexicographical order used to compute linear indices

temp=sort(newLexOrder);
if(~isequal(temp,1:A.order()))
    error('For N degree tensor, lexicographical order must be a permutation of the numbers within 1:N');
end

if(~isequal(A.lexOrder,newLexOrder))
    A.indices=A.indices-1;
    SparseChangeLinIdxMex(A.indices,uint8(A.lexOrder),uint8(newLexOrder),int64(A.dims));
    A.lexOrder=newLexOrder;
    A.indices=A.indices+1;
    A.isSorted=false;
end

end

