function is_eq = eq(A,B)
%EQ Summary of this function goes here
%   Detailed explanation goes here
func=@isequal;
if(isa(A,'SparseMIA'))
    is_eq= compare(A,B,func);
else
    is_eq= compare(B,A,func);

end

