function is_eq = eq(A,B)
%EQ Summary of this function goes here
%   Detailed explanation goes here
func=@isequal;
is_eq= compare(A,B,func);

end

