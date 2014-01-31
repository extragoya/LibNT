function C = intersect(A,B,fhandle)
% INTERSECT (intersect sub-operation for sparse lattices)

%Originally programmed by Dileepan Joseph, updated by Adam Harrison
%UofA Electronic Imaging Lab 2010

% executes binary operations analagous to the intersect function (and)

nnzA = A.nnz;
nnzB = B.nnz;
c_inds=zeros(nnzA);
c_vals=zeros(nnzA); %worst case, all nnz values intersect
a_ctr = 1;
b_ctr = 1;
c_ctr = 1; %combined index

while a_ctr <= nnzA && b_ctr <= nnzB
    if A.inds(a_ctr) < B.inds(b_ctr)
        
        a_ctr = a_ctr+1;
        
    elseif A.inds(a_ctr) > B.inds(b_ctr)
        
        b_ctr = b_ctr+1;
        
    else
        c_inds(c_ctr)=A.inds(a_ctr);        
        c_vals(c_ctr)=fhandle(A.vals(a_ctr),B.vals(b_ctr));
        a_ctr = a_ctr+1;
        b_ctr = b_ctr+1;
        c_ctr=c_ctr+1;
    end
    
end

c_inds=c_inds(1:c_ctr-1);   %remove extra space
c_vals=c_vals(1:c_ctr-1);

C = SparseLattice(c_inds,c_vals,A.m,A.n,A.p);
