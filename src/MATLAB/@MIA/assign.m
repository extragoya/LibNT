function obj = assign(obj,otherMIA,assign_order)
%assigns from otherMIA to obj, based on the index_order, e.g., [2 1 3]
%means that otherMIA is assigned to obj as obj(ijk)=otherMIA(jik)

if(isa(otherMIA,'MIA'))
    obj.data=permute(otherMIA.data,assign_order);
    
else %TODO as SparseMIA (above MIA)
    error('Incompatible object used in assignment');
end

end

