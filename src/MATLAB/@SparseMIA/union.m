function C = union(A,B,fhandle)
% UNION (union sub-operation for sparse mias)


% executes binary operations analagous to the union function (add,
% subtract, or)



nnzA = size(A.data,1);
nnzB = size(B.data,1);
c_inds=zeros(nnzA+nnzB,length(A.partition));  %in worst case, no intersection, so C has twice as many nnz
c_vals=zeros(nnzA+nnzB,1);

a_ctr = 1;
b_ctr = 1;
c_ctr = 1; %combined index

exponent=size(A.indices,2)-1:-1:0;
mult=ones(1,size(A.indices,2))*2;
mult=mult.^exponent;

while a_ctr <= nnzA && b_ctr <= nnzB
    comparison=(A.indices(a_ctr,:)-B.indices(b_ctr,:));
    comparison(comparison<0)=-1;
    comparison(comparison>0)=1;
    comparison=sum(comparison.*mult);
    if comparison < 0        
        c_inds(c_ctr,:)=A.indices(a_ctr,:);
        c_vals(c_ctr)=fhandle(A.data(a_ctr),0);
        a_ctr = a_ctr+1;
        
    elseif comparison > 0        
        c_inds(c_ctr,:)=B.indices(b_ctr,:);
        c_vals(c_ctr)=fhandle(0,B.data(b_ctr));
        b_ctr = b_ctr+1;
        
    else        
        c_inds(c_ctr,:)=A.indices(a_ctr,:);
        c_vals(c_ctr)=fhandle(A.data(a_ctr),B.data(b_ctr));       
        a_ctr = a_ctr+1;
        b_ctr = b_ctr+1;
       
    end
    c_ctr=c_ctr+1; 
end

if a_ctr<=nnzA
    c_inds(c_ctr:c_ctr+(nnzA-a_ctr),:)=A.indices(a_ctr:end,:);
    c_vals(c_ctr:c_ctr+(nnzA-a_ctr))=fhandle(A.data(a_ctr:end),0);
    c_ctr=c_ctr+(nnzA-a_ctr);
elseif b_ctr<=nnzB
    c_inds(c_ctr:c_ctr+(nnzB-b_ctr),:)=B.indices(b_ctr:end,:);
    c_vals(c_ctr:c_ctr+(nnzB-b_ctr))=fhandle(0,B.data(b_ctr:end));
    c_ctr=c_ctr+(nnzB-b_ctr);
else
    c_ctr=c_ctr-1;
end

c_inds=c_inds(1:c_ctr,:);   %remove extra space
c_vals=c_vals(1:c_ctr);

C = SparseMIA(c_vals,c_inds,A.dims,A.partition);