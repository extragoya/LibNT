function C= sparse_times(A,B)
%Multiplying a delta lattice, i.e., a lattice constructed from a Dirac
%delta NT, with a sparse lattice. Instead of explicitly multiplying we can
%directly construct the product


  


if(isa(A,'DeltaLattice'))
    delta_dim=A.delta_dim;
    if(A.row_idx==2) %if we're performing a pure  outer product
       c_inds=B.inds;
       c_vals=B.vals;
       delta_inds=int64(0:delta_dim-1);
       delta_inds=delta_inds';
       delta_inds=delta_inds+delta_inds*delta_dim;
       delta_inds=repmat(delta_inds,B.nnz,1);
       c_inds=c_inds';
       c_inds=repmat(c_inds,delta_dim,1);
       c_inds=c_inds-1;
       c_inds=c_inds*delta_dim^2;
       c_inds=c_inds(:)+delta_inds;
       c_vals=c_vals';
       c_vals=repmat(c_vals,delta_dim,1); 
       c_vals=c_vals(:);
       c_inds=c_inds+1;
       C=SparseLattice(c_inds,c_vals,delta_dim^2,B.n,B.p);


    end
    
else
    delta_dim=B.delta_dim;
    if(B.col_idx==2) %if we're performing a pure  outer product
       c_inds=A.inds;
       c_vals=A.vals;
       delta_inds=int64(0:delta_dim-1);       
       delta_inds=delta_inds+delta_inds*delta_dim; %create indices from delta lattice
       delta_inds=repmat(delta_inds,A.nnz,1); %replicate them - each index is matched to every index in A
       delta_inds=delta_inds(:);
       delta_inds=delta_inds*A.m; %since the delta indices will have the most significance, multiply them by B.m
       
       c_inds=repmat(c_inds,delta_dim,1);
       c_inds=c_inds-1;       
       c_inds=c_inds+delta_inds;       
       c_vals=repmat(c_vals,delta_dim,1);        
       c_inds=c_inds+1;
       C=SparseLattice(c_inds,c_vals,A.m,delta_dim^2,B.p);

    elseif(B.col_idx==1 && B.depth_idx==1) %performing an outer and element-wise product
        c_inds=A.inds-1;
        c_depth_inds=idivide(c_inds,A.m);
        c_inds=c_inds+c_depth_inds*A.m*A.p;
        c_inds=c_inds+1;
        C=SparseLattice(c_inds,A.vals,A.m,delta_dim,A.p);
    end
end





end

