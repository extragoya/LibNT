function is_eq = eq(A,B)
%EQ Summary of this function goes here
%   Detailed explanation goes here
if ~isa(A,'DenseNT') || ~isa(B,'DenseNT')
    error('SparseNTs can only be multipled with NT classes');
end
func=@isequal;
if(isa(A,'SparseNT'))
    is_eq= compare(A,B,func);
else
    is_eq= compare(B,A,func);

end

