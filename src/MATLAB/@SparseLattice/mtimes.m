function C=mtimes(A,B)
% mtimes (multiply two sparse lattices)

%Function performs lattice multiplication by creating a sparse matrix for
%each tab in the two lattices. The sparse matrices are created by finding
%the unique rows and columsn of the current tabs of A and B, respectively.
%As well, the intersection of unique columns and rows of A and B are also
%found, respectively. These indices are used to map the current tabs into
%sparse matrices that have NO non-zero columns and rows. As MATLAB uses
%compressed column format, this allows one to mulitply sparse lattices with
%huge numbers of columns or rows but relatively small numbers of nonzeros.

%Adam Harrison
%UofA Electronic Imaging Lab 2011


if A.n~=B.m
    error('Number of columns in A must equal number of rows in B')
elseif A.p~=B.p
    error('Depth of both lattices must be identical')
end
m=A.m;
n=B.n;
p=A.p;
dimen=m*n*p;
if isa(B,'SparseLattice') && isa(A,'SparseLattice')
    A.inds=A.inds-1;
    B.inds=B.inds-1;
    [c_data, c_inds]=SparseLatticeMultMex(A.vals,int64(A.inds),[A.m A.n A.p],B.vals,int64(B.inds),[B.m B.n B.p]);
    c_inds=c_inds+1;
    C=SparseLattice(c_inds,c_data,m,n,p);
    A.inds=A.inds+1;
    B.inds=B.inds+1;
else
    if isa(A,'SparseLattice')
        A.inds=A.inds-1;
        [c_data, c_inds]=SparseDenseLatticeMultMex(A.vals,int64(A.inds),[A.m A.n A.p],B.vals,0);
        A.inds=A.inds+1;
    else
        B.inds=B.inds-1;
        [c_data, c_inds]=SparseDenseLatticeMultMex(B.vals,int64(B.inds),[B.m B.n B.p],A.vals,1);
        B.inds=B.inds+1;
    end 
    c_inds=c_inds+1;
    if length(c_data)/dimen <0.5
        C=SparseLattice(c_inds,c_data,m,n,p);
    else
        c_vals=zeros(m,n,p);
        c_vals(c_inds)=c_data;
        C=Lattice(c_vals);
    end
    
    
end









