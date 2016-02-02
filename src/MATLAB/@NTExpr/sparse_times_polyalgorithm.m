function CExpr=sparse_times_polyalgorithm(A,B,a_inner_idx,a_outer_idx,a_inter_idx,b_outer_idx,b_inner_idx, b_inter_idx)
%This implements the sparse multiplication poly-algorithm which dispatches
%to different multiplication routines depending on the hyper-sparsity of
%the operands
A_nt=A.m_mia;
B_nt=B.m_mia;
a_dims=A_nt.dims();

a_inter_dims=a_dims(a_inter_idx);
a_inner_dims=a_dims(a_inner_idx);
a_inner_size=prod(a_inner_dims);
a_col_sparsity=(a_inner_size/A_nt.nnz) >3;

a_outer_dims=a_dims(a_outer_idx);
a_outer_size=prod(a_outer_dims);
a_row_sparsity=(a_outer_size/A_nt.nnz) >3;

a_sparse=~(a_col_sparsity || a_row_sparsity);
a_hyper_sparse=a_col_sparsity && a_row_sparsity;
a_row_sparsity=a_row_sparsity&~a_hyper_sparse;
a_col_sparsity=a_col_sparsity&~a_hyper_sparse;

b_dims=B_nt.dims();
b_inner_dims=b_dims(b_inner_idx);
b_inner_size=prod(b_inner_dims);
b_row_sparsity=(b_inner_size/B_nt.nnz) >3;

b_outer_dims=b_dims(b_outer_idx);
b_outer_size=prod(b_outer_dims);
b_col_sparsity=(b_outer_size/B_nt.nnz) >3;

b_sparse=~(b_col_sparsity || b_row_sparsity);
b_hyper_sparse=b_col_sparsity && b_row_sparsity;
b_row_sparsity=b_row_sparsity&~b_hyper_sparse;
b_col_sparsity=b_col_sparsity&~b_hyper_sparse;

if a_sparse
    if b_sparse %both are normal sparse, so just do CSC or CSR
        if(b_outer_size>a_outer_size)
            transposed=0;
            algorithm=1;
            
        else
            transposed=1;
            algorithm=1;
        end
    else    %if A is sparse, but B has some hyper-sparsity, do CSC
       transposed=0;
        algorithm=1; 
    end
elseif b_sparse %in this case B is sparse, but A has some hyper-sparsity, so do CSR
    transposed=1;
    algorithm=1;
elseif a_col_sparsity
    if b_row_sparsity %if A is col-sparse and B is row-sparse, we do DCSC or DCSR
        if(b_outer_size>a_outer_size)
            transposed=0;
            algorithm=2;
            
        else
            transposed=1;
            algorithm=2;
        end
    else %A is col-sparse and B is col-sparse or hyper-sparse. so do DCSC
        transposed=0;
        algorithm=2;
    end
elseif b_row_sparsity %B is row sparse and A is row- or hyper-sparse, so do DCSR
    transposed=1;
    algorithm=2;
elseif a_row_sparsity %A is row sparse, and B is col-sparse or hyper-sparse
    if b_col_sparsity %if B is col-sparse do CSCNA or CSRNA
        if(b_outer_size>a_outer_size)
            transposed=0;
            algorithm=3;
            
        else
            transposed=1;
            algorithm=3;
        end
    else %otherwise if B is hyper-sparse do CSCNA
        transposed=0;
        algorithm=3;
    end
elseif b_col_sparsity %B is column-sparse and A is hyper-sparse, so do CSRNA
   transposed=1;
   algorithm=3; 
else %A and B are both hyper-sparse, so we do outer-product algorithm
   transposed=2;
   algorithm=4; 
end

if transposed==0
    
    A_lat=A_nt.toLattice(a_outer_idx, a_inner_idx, a_inter_idx);
    B_lat=B_nt.toLattice(b_inner_idx, b_outer_idx, b_inter_idx);
    %perform lattice product
    C_lat=specialized_mult(A_lat,B_lat,algorithm);
    
    %convert result back to an MIA
    c_indices=A.m_indices(a_outer_idx);
    c_indices=[c_indices B.m_indices(b_outer_idx)];
    c_indices=[c_indices A.m_indices(a_inter_idx)];

    C=C_lat.toMIA(a_outer_dims, b_outer_dims, a_inter_dims);
    
    
elseif transposed==1
    
    A_lat=A_nt.toLattice(a_inner_idx,a_outer_idx,a_inter_idx);
    B_lat=B_nt.toLattice(b_outer_idx,b_inner_idx, b_inter_idx);
    %perform lattice product
    
    C_lat=specialized_mult(B_lat,A_lat,algorithm);
    
    
    
    %convert result back to an MIA
    c_indices=A.m_indices(a_outer_idx);
    c_indices=[c_indices B.m_indices(b_outer_idx)];
    c_indices=[c_indices A.m_indices(a_inter_idx)];

    C=C_lat.toMIA(b_outer_dims,a_outer_dims, a_inter_dims,true);
    
else %outer product version
    A_lat=A_nt.toLattice(a_outer_idx, a_inner_idx, a_inter_idx);
    B_lat=B_nt.toLattice(b_outer_idx,b_inner_idx, b_inter_idx); %a transposed version of B, specialized_mult will turn it to a Rowmajor lattice
     %perform lattice product
    C_lat=specialized_mult(A_lat,B_lat,algorithm);
    
    %convert result back to an MIA
    c_indices=A.m_indices(a_outer_idx);
    c_indices=[c_indices B.m_indices(b_outer_idx)];
    c_indices=[c_indices A.m_indices(a_inter_idx)];

    C=C_lat.toMIA(a_outer_dims, b_outer_dims, a_inter_dims);
end

%create a MIAExpr using resulting MIA and indices
CExpr=NTExpr(C,c_indices,length(a_inter_idx));

end

