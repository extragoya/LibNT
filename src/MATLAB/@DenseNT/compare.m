function is_equal=compare(A,B,func)
%compare two NTs with each other, func can be equals or a fuzzy compare
is_equal=true;
if(isa(B,'DenseNT'))
    if (~isequal(A.dims,B.dims))
        is_equal=false;
    end    
else 
    error('Incompatible object used in comparison'); 
end

is_equal=is_equal&& func(A.data,B.data);

end

