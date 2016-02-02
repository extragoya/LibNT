function C=concat(A,B)
%concats two sparseNTs together to create order+1 NT. If old dims were
%[DIMS], new dims are [2 DIMS]

if ~isa(A,'SparseNT') || ~isa(B,'SparseNT')
    error('Only concat two sparseNTs together');
end

if A.dims()~=B.dims()
    error('Can only concat two SparseNTs together of the same dimensionality and degree')
end
lexOrder=1:A.order();

A=A.permute(lexOrder);
B=B.permute(lexOrder);
b_inds=B.indices;
b_inds=b_inds+A.dimensionality();
C=SparseNT([A.data;B.data],[A.indices;b_inds],[2 A.dims()],[lexOrder+1 1],true);

end

