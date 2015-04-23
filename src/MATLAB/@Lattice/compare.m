function is_eq = compare(A,B,func)
%EQ Summary of this function goes here
%   Detailed explanation goes here
is_eq= false;
if ~isa(B,'Lattice')
    return;
end

if(A.m~=B.m || A.m~=B.m || A.m~=B.m)
    return;
end

is_eq=func(A.vals,B.vals);


end

