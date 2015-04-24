function obj = assign(obj,otherMIA,assign_order)
%assigns from otherMIA to obj, based on the index_order, e.g., [2 1 3]
%means that otherMIA is assigned to obj as obj(ijk)=otherMIA(jik)

if(isa(otherMIA,'SparseMIA') && isa(obj,'SparseMIA'))
    obj.data=otherMIA.data;
    obj.indices=otherMIA.indices; 
    
    obj.mDims=otherMIA.dims;
    obj.mDims=obj.mDims(assign_order);
    reverse_order=zeros(1,length(assign_order));
    reverse_order(assign_order)=1:length(assign_order);
    %obj.lexOrder=otherMIA.lexOrder(reverse_order);
    obj.lexOrder=reverse_order(otherMIA.lexOrder);
    if(~otherMIA.isSorted)
        obj.isSorted=false;
    elseif(~isequal(assign_order',1:length(assign_order)))
        obj.isSorted=false;
    else
        obj.isSorted=true;
        
    end
    
else %TODO as MIA (above SparseMIA)
    error('Incompatible object used in assignment');
end

end

