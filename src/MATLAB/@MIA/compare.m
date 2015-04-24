function is_equal=compare(A,B,func)
%assigns from otherMIA to obj, based on the index_order, e.g., [2 1 3]
%means that otherMIA is assigned to obj as obj(ijk)=otherMIA(jik)
is_equal=true;
if(isa(B,'MIA'))
    if (~isequal(A.dims,B.dims))
        is_equal=false;
    end    
else 
    error('Incompatible object used in comparison'); 
end

is_equal=is_equal&& func(A.data,B.data);

end

