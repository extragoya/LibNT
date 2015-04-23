function C= do_merge(A,B,permute_idx,op)
%perform addition based on the permutation idx of B
if ~isa(A,'SparseMIA') || ~isa(B,'SparseMIA')
    error('do_merge can only be called with classes of SparseMIA');
end

permuteA=false;

if(A.isSorted || B.isSorted)
    
    if(~A.isSorted)
        permuteA=true;    
    elseif (B.isSorted && A.nnz<B.nnz)
        permuteA=true;         
    end
else
    if(A.nnz>B.nnz)
        A=A.sort();        
    else
        B=B.sort();
        permuteA=true;
    end
end

if(permuteA)
    reverse_order=zeros(1,length(permute_idx));
    reverse_order(permute_idx)=1:length(permute_idx);
    newLexOrder=B.lexOrder;
    newLexOrder=reverse_order(newLexOrder);
    A=A.permute(newLexOrder);    
else  
    newLexOrder=A.lexOrder;
    newLexOrder=permute_idx(newLexOrder)';
    B=B.permute(newLexOrder);    
end
    
[c_data, c_inds]=SparseMergeMex(A.data,A.indices,B.data,B.indices,op);

C=SparseMIA(c_data,c_inds,A.dims,A.lexOrder,true);

end

