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
    [c_data c_inds]=sparselattice_mul(A.vals,A.inds,[A.m A.n A.p],B.vals,B.inds,[B.m B.n B.p]);
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

% 
% a_end=1;
% b_end=1;
% %used to store values of C and indices
% final_row_inds=[];
% final_col_inds=[];
% final_depth_inds=[];
% final_vals=[];
% 
% for current_depth=1:p %treat each depth separately
%     need_update=false;
%     %find entries belonging to current depth
%     a_begin=a_end;
%     while A.depth(a_end)==current_depth
%         a_end=a_end+1;
%     end
%     
%     if isa(B,'SparseLattice')
%         b_begin=b_end;
%         while B.depth(b_end)==current_depth
%             b_end=b_end+1;
%         end
%         %only need to multiply current tab if both lattices have nonzeros in
%         %current depth
%         if b_end~=b_begin && a_end~=a_begin
%             need_update=true;
%             a_col=A.col(a_begin:a_end-1);
%             b_col=B.col(b_begin:b_end-1);
%             a_row=A.row(a_begin:a_end-1);
%             b_row=B.row(b_begin:b_end-1);
%             
%             %find unique columns and rows of A and B
%             [a_col_unique,~,a_col_map] = unique(a_col);
%             
%             [b_row_unique,~,b_row_map] = unique(b_row);
%             
%             %find their intersection
%             [~, ia, ib] = intersect(a_col_unique, b_row_unique);
%             length_inner=length(ia);
%             
%             %find which entries in A belong to the intersection
%             tf=false(length(a_col_unique),1);
%             tf(ia)=true;
%             idx=tf(a_col_map);
%             %extract only intersected values
%             a_vals=A.vals(a_begin:a_end-1);
%             a_vals=a_vals(idx);
%             a_row=a_row(idx);
%             a_col_map=a_col_map(idx);
%             [~,~,a_col_map] = unique(a_col_map); %update A column mapping post intersection
%             clear ia idx
%             
%             %do the same to B
%             tf=false(length(b_row_unique),1);
%             tf(ib)=true;
%             idx=tf(b_row_map);
%             b_vals=B.vals(b_begin:b_end-1);
%             b_vals=b_vals(idx);
%             b_col=b_col(idx);
%             b_row_map=b_row_map(idx);
%             [~,~,b_row_map] = unique(b_row_map);
%             clear ib idx
%             
%             %find unique rows and columns of A and B
%             [a_row_unique,~,a_row_map] = unique(a_row);
%             
%             [b_col_unique,~,b_col_map] = unique(b_col);
%             
%             %create sparse matrices A and B based on the mappings
%             A_matrix=sparse(a_row_map,a_col_map,a_vals,length(a_row_unique),length_inner);
%             B_matrix=sparse(b_row_map,b_col_map,b_vals,length_inner,length(b_col_unique));
%             clear a_vals b_vals a_row_map a_col_map b_row_map b_col_map
%             %compute the product
%             C_matrix=A_matrix*B_matrix;
%             clear A_matrix B_matrix
%             %extract values and compressed indices of C matrix
%             [ii jj s]=find(C_matrix);
%             clear C_matrix
%             %map them back to their uncompressed indices
%             row_inds=a_row_unique(ii);
%             col_inds=b_col_unique(jj);
%             depth_inds=ones(length(row_inds),1)*current_depth;
%             
%             
%         end
%     else
%         if a_end~=a_begin
%             need_update=true;
%             a_vals=A.vals(a_begin:a_end-1);
%             a_col=A.col(a_begin:a_end-1);
%             a_row=A.row(a_begin:a_end-1);
%             
%             %find unique columns of A
%             [a_col_unique,~,a_col_map] = unique(a_col);
%             %only extract those columns from B
%             B_matrix=B(a_col_unique,:,current_depth);
%             
%             %find unique rows of A
%             [a_row_unique,~,a_row_map] = unique(a_row);
%             %create sparse matrices A and B based on the mappings
%             A_matrix=sparse(a_row_map,a_col_map,a_vals,length(a_row_unique),length(a_col_unique));
%             
%             clear a_vals a_row_map a_col_map
%             %compute the product
%             C_matrix=A_matrix*B_matrix;
%             clear A_matrix B_matrix
%             %extract values and compressed indices of C matrix
%             [ii jj s]=find(C_matrix);
%             clear C_matrix
%             %map them back to their uncompressed indices
%             row_inds=a_row_unique(ii);
%             col_inds=jj;
%             depth_inds=ones(length(row_inds),1)*current_depth;
%             
%         end
%         
%         
%         
%         
%     end
%     
%     %update the list of values and indices for C lattice
%     if need_update
%         if isempty(final_row_inds)
%             final_row_inds=row_inds;
%             final_col_inds=col_inds;
%             final_depth_inds=depth_inds;
%             final_vals=s;
%         else
%             final_length=length(final_row_inds);
%             current_length=length(row_inds);
%             temp=zeros(final_length+current_length,1);
%             
%             temp(1:final_length)=final_row_inds;
%             temp(final_length+1:end)=row_inds;
%             final_row_inds=temp;
%             
%             temp(1:final_length)=final_col_inds;
%             temp(final_length+1:end)=col_inds;
%             final_col_inds=temp;
%             
%             temp(1:final_length)=final_depth_inds;
%             temp(final_length+1:end)=depth_inds;
%             final_depth_inds=temp;
%             
%             temp(1:final_length)=final_vals;
%             temp(final_length+1:end)=s;
%             final_vals=temp;
%         end
%     end
%     
% end
%create the sparse lattice from the list of values







