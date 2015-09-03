function C = mrdivide(A,B)
%Scalar multiplication


if(isscalar(B) && isa(A,'DenseNT'))
    C=A;
    C.data=C.data/B;
else
   error('NTs can only be divided with scalars using ./ operator. If you are trying to perform NT multiplication use index notation, e.g., A(''ij'')*B(''jk'')');
end







end

