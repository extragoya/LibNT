function C=mtimes(A,B)
% mtimes (multiply two sparse lattices)

%Adam Harrison
%UofA Electronic Imaging Lab 2010

%Naively, one would compute C=A*B. In such a configuration, A would be
%transposed first in order to travel down its rows (as SparseLattice
%indices are column major). But in doing so, C would be ordered in a
%row-major sense, requiring a subsequent sort. If C has many more nnzs than
%A and B, it's best to avoid a potentially costly sort.

%So function actually computes C=(B'*A')', where the transpose operation
%occurs in-line. By swapping row and columns of (B'*A') values to fill in C
%the function performs a transpose and fills in C in column-major ordering.



if A.n~=B.m
    error('Number of columns in A must equal number of rows in B')
elseif A.p~=B.p
    error('Depth of both lattices must be identical')
end


b_n=B.n;
B=B'; %we want B' not B in order to get row-major ordering for B






a_m=A.m;
a_n=A.n;
a_p=A.p;
nnzA=A.nnz;


nnzB=B.nnz;

C=SparseLattice(a_m,b_n,a_p);
N0=nnzB+nnzA;

a_ctr=1;
b_ctr=1;


f=[];





for current_depth=1:a_p
    
    
    while A.depth(a_ctr)==current_depth && B.depth(b_ctr)==current_depth
        b_row=B.col(b_ctr); %the row of B is the col of B'
        a_col=A.col(a_ctr);
        if a_col<b_row
            a_ctr=a_ctr+1;
        elseif b_row<a_col
            b_ctr=b_ctr+1;
        else
                      
            N_temp=N0;
            temp_f=zeros(N0,1);
            temp_f_inds=temp_f;
            temp_f_ctr=0;
            b_ctr2=b_ctr;
            while b_ctr2<=nnzB && B.col(b_ctr2)==b_row
                a_ctr2=a_ctr;
                while a_ctr2<=nnzA &&  A.col(a_ctr2)==a_col
                    temp_f_ctr=temp_f_ctr+1;
                    temp_f(temp_f_ctr)=A.vals(a_ctr2)*B.vals(b_ctr2);
                    temp_f_inds(temp_f_ctr)=sub2ind(C,A.row(a_ctr2),B.row(b_ctr2),current_depth);
                    
                    if temp_f_ctr>N_temp;
                        temp_f=[temp_f;zeros(N_temp,1)];
                        N_temp=N_temp*2;
                    end
                    
                    
                    a_ctr2=a_ctr2+1;
                end
                b_ctr2=b_ctr2+1;
            end
            
            %if this is the first list of numbers to be generated, store them
            if isempty(f);
                f=temp_f(1:temp_f_ctr);
                f_inds=temp_f_inds(1:temp_f_ctr);                
            else %otherwise, we need sort the current list with the running list in 'f'
                
                [f f_inds]=SparseLattice.merge_list(f,f_inds,temp_f(1:temp_f_ctr),temp_f_inds(1:temp_f_ctr),a_m,b_n,a_p);
                
            end
            
            
            
            a_ctr=a_ctr2;
            b_ctr=b_ctr2;
            
        end
        
        
    end
    %move counters forward until next lattice tab
    while A.depth(a_ctr)==current_depth
        a_ctr=a_ctr+1;
    end
    while B.depth(b_ctr)==current_depth
        b_ctr=b_ctr+1;
    end
    
    
end


C.inds=f_inds;
C.vals=f;
C.nnz=length(f_inds);

