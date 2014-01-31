function C = union(A,B,fhandle)
% UNION (union sub-operation for sparse lattices)

%Originally programmed by Dileepan Joseph, updated by Adam Harrison
%UofA Electronic Imaging Lab 2010

% executes binary operations analagous to the union function (add,
% subtract, or)



nnzA = A.nnz;
nnzB = B.nnz;
c_inds=zeros(nnzA+nnzB,1);  %in worst case, no intersection, so C has twice as many nnz
c_vals=zeros(nnzA+nnzB,1);

a_ctr = 1;
b_ctr = 1;
c_ctr = 1; %combined index



while a_ctr <= nnzA && b_ctr <= nnzB
    
    if A.inds(a_ctr) < B.inds(b_ctr)        
        c_inds(c_ctr)=A.inds(a_ctr);
        c_vals(c_ctr)=fhandle(A.vals(a_ctr),0);
        a_ctr = a_ctr+1;
        
    elseif A.inds(a_ctr) > B.inds(b_ctr)        
        c_inds(c_ctr)=B.inds(b_ctr);
        c_vals(c_ctr)=fhandle(0,B.vals(b_ctr));
        b_ctr = b_ctr+1;
        
    else        
        c_inds(c_ctr)=A.inds(a_ctr);
        c_vals(c_ctr)=fhandle(A.vals(a_ctr),B.vals(b_ctr));       
        a_ctr = a_ctr+1;
        b_ctr = b_ctr+1;
       
    end
    c_ctr=c_ctr+1; 
end

if a_ctr<=nnzA
    c_inds(c_ctr:c_ctr+(nnzA-a_ctr))=A.inds(a_ctr:end);
    c_vals(c_ctr:c_ctr+(nnzA-a_ctr))=fhandle(A.vals(a_ctr:end),0);
    c_ctr=c_ctr+(nnzA-a_ctr);
elseif b_ctr<=nnzB
    c_inds(c_ctr:c_ctr+(nnzB-b_ctr))=B.inds(b_ctr:end);
    c_vals(c_ctr:c_ctr+(nnzB-b_ctr))=fhandle(0,B.vals(b_ctr:end));
    c_ctr=c_ctr+(nnzB-b_ctr);
else
    c_ctr=c_ctr-1;
end

c_inds=c_inds(1:c_ctr);   %remove extra space
c_vals=c_vals(1:c_ctr);

C = SparseLattice(c_inds,c_vals,A.m,A.n,A.p);
