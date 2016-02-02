function obj = assign(obj,otherNT,assign_order)
%assigns from otherMIA to obj, based on the index_order, e.g., [2 1 3]
%means that otherMIA is assigned to obj as obj(ijk)=otherMIA(jik)

if(isa(otherNT,'SparseNT') && isa(obj,'SparseNT'))
    obj.data=otherNT.data;
    obj.indices=otherNT.indices; 
    
    obj.mDims=otherNT.dims;
    obj.mDims=obj.mDims(assign_order);
    %creating the newLexOrder requires using the inverse shuffle of assign
    %order
    reverse_order=zeros(1,length(assign_order));
    reverse_order(assign_order)=1:length(assign_order);
    
    obj.lexOrder=reverse_order(otherNT.lexOrder);
    if(~otherNT.isSorted)
        obj.isSorted=false;
    
    else
        obj.isSorted=true;
        
    end
    
else %TODO as DenseNT (above SparseNT)
    error('Incompatible object used in assignment');
end

end

