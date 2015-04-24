function is_equal=compare(A,B,func)
%assumes A is SparseMIA. Accepts B to be MIA or SparseMIA
is_equal=true;
if(isa(B,'MIA'))
    if (~isequal(A.dims,B.dims))
        is_equal=false;
    end    
else 
    error('Incompatible object used in comparison'); 
end
if (is_equal)
    if(isa(B,'SparseMIA'))
        B.changeLexOrder(A.lexOrder);
        is_equal=isequal(A.indices,B.indices);
        is_equal=is_equal && func(A.data,B.data);
    
    elseif (isa(B,'MIA')) %TODO as MIA (above SparseMIA)
        A.changeLexOrder(1:length(A.dims));
        B_inds=int64(find(B.data));
        B_vals=B.data(B_inds);
        B_vals=B_vals(:);
        if(isempty(B_vals) && isempty(A.vals))
            is_equal=true;
        elseif(length(A.data)~=length(B_vals))
        
            is_equal=false;
        elseif(isequal(A.indices,B_inds) && func(A.data,B_vals))
        
            is_equal= true;
        end

    end
end
end

