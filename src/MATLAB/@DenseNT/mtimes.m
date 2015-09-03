function C = mtimes(A,B)
%Scalar multiplication

if(isscalar(A))
    C=B;
    C.data=C.data*A;
elseif(isscalar(B))
    C=A;
    C.data=C.data*B;
else
   error('NTs can only be multiplied with scalars using * operator. If you are trying to perform NT multiplication use index notation, e.g., A(''ij'')*B(''jk'')');
end







end

