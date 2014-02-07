function is_eq = compare(A,B,func)
%EQ Summary of this function goes here
%   Detailed explanation goes here
is_eq= false;
if isa(B,'Lattice')
    if(A.m~=B.m || A.m~=B.m || A.m~=B.m)
        return;
    end
end

if isa(B,'SparseLattice')
    if(isequal(A.inds,B.inds) && func(A.vals,B.vals))
        
        is_eq= true;
    end
elseif isa(B,'Lattice')
    
    B_inds=int64(find(B.vals));
    B_vals=B.vals(B_inds);
    B_vals=B_vals(:);
    if(length(A.vals)~=length(B_vals))
        
        is_eq=false;
    elseif(isequal(A.inds,B_inds) && func(A.vals,B_vals))
        
        is_eq= true;
    end
    
end


end

