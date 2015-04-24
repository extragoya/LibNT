function C = mrdivide(A,B)
%Scalar multiplication


if(isscalar(B) && isa(A,'MIA'))
    C=A;
    C.data=C.data/B;
else
   error('MIAs can only be multiplied with scalars using * operator. If you are trying to perform MIA multiplication use index notation, e.g., A(''ij'')*B(''jk'')');
end







end

