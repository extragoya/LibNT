function A=transpose(A)
% transpose for sparse lattices

%Adam Harrison
%UofA Electronic Imaging Lab 2010

%executes transpose of each matrix within lattice


[i j k]=ind2sub(A);

a_n=A.n;
A.n=A.m;
A.m=a_n;

inds=sub2ind(A,j,i,k); %flip the rows and column indices around

[inds idx]=sort(inds); %sort them, and return transposed version of A
A.inds=inds;
A.vals=A.vals(idx);