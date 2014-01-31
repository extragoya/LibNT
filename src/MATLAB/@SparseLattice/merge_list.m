function [c c_inds]=merge_list(a,a_inds,b,b_inds,m,n,p)

if length(a)~=length(a_inds)
    error('a and a_inds must have identical lengths')
end
if length(b)~=length(b_inds)
    error('b and b_inds must have identical lengths')
end

A=SparseLattice(a_inds,a,m,n,p);
B=SparseLattice(b_inds,b,m,n,p);
C=A+B;
c=C.vals;
c_inds=C.inds;

