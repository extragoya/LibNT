function obj = assign(obj,otherNT,assign_order)
%assigns from otherNT to obj, based on the index_order, e.g., [2 1 3]
%means that otherNT is assigned to obj as obj(ijk)=otherNT(jik)
if (isa(otherNT,'SparseNT'))
    obj=SparseNT;
    obj=obj.assign(otherNT,assign_order);
elseif(isa(otherNT,'DenseNT') )
    if otherNT.order>1
        obj.data=permute(otherNT.data,assign_order);
    else
        obj.data=otherNT.data;
    end
    
else %TODO as SparseNT (above NT)
    error('Incompatible object used in assignment');
end

end

