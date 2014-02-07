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
if isa(B,'SparseLattice')
    A.inds=A.inds-1;
    B.inds=B.inds-1;
    [c_data, c_inds]=SparseLatticeMultMex(A.vals,A.inds,[A.m A.n A.p],B.vals,B.inds,[B.m B.n B.p]);
    c_inds=c_inds+1;
    C=SparseLattice(c_inds,c_data,m,n,p);
    A.inds=A.inds+1;
    B.inds=B.inds+1;
else
    a_end=1;
    final_row_inds=[];
    final_col_inds=[];
    final_depth_inds=[];
    final_vals=[];
    for current_depth=1:p %treat each depth separately
        a_begin=a_end;
        while a_end<length(A.vals) && A.tab(A.inds(a_end))==current_depth
            a_end=a_end+1;
        end
        if a_end~=a_begin
            
            a_vals=A.vals(a_begin:a_end-1);
            a_col=A.col(A.inds(a_begin:a_end-1));
            a_row=A.row(A.inds(a_begin:a_end-1));
            
            
            B_matrix=B(:,:,current_depth);
            
            
            A_matrix=sparse(double(a_row),double(a_col),a_vals,A.m,A.n);
            
            C_matrix=A_matrix*B_matrix;
            clear A_matrix B_matrix
            
            [ii jj s]=find(C_matrix);
            clear C_matrix
            %map them back to their uncompressed indices
            
            depth_inds=ones(length(ii),1)*current_depth;
            
            
            if isempty(final_row_inds)
                final_row_inds=ii;
                final_col_inds=jj;
                final_depth_inds=depth_inds;
                final_vals=s;
            else
                final_length=length(final_row_inds);
                current_length=length(ii);
                temp=zeros(final_length+current_length,1);
                
                temp(1:final_length)=final_row_inds;
                temp(final_length+1:end)=ii;
                final_row_inds=temp;
                
                temp(1:final_length)=final_col_inds;
                temp(final_length+1:end)=jj;
                final_col_inds=temp;
                
                temp(1:final_length)=final_depth_inds;
                temp(final_length+1:end)=depth_inds;
                final_depth_inds=temp;
                
                temp(1:final_length)=final_vals;
                temp(final_length+1:end)=s;
                final_vals=temp;
            end
            
            
        end
        
        
    end
    C=SparseLattice(final_row_inds,final_col_inds,final_depth_inds,final_vals,m,n,p);
    
end









