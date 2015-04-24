function C=mtimes(A,B)

if ~isa(A,'MIA') || ~isa(B,'MIA')
    error('DeltaMIAs can only be multipled with classes or subclasses of MIA');
end

error_check_mult(A,B); %make sure indices match up properly for multiplication

flipped=false;
if ~isa(A,'DeltaMIA')
    flipped=true;
    temp=B;
    B=A;
    A=temp;
    clear temp
end


dims=size(B);


if ~isa(B,'SparseMIA')
    B=SparseMIA(B);

end
a_inter=A.inter_idx;
b_inter=B.inter_idx;
b_outer=B.outer_idx;
b_inner=B.inner_idx;
if ~isempty(a_inter)
    if ~isempty(A.inner_idx) %inter/inner product
        idx1=B.extract_indices(b_inter);
        idx2=B.extract_indices(b_inner);
        idx1=idx2-idx1;
        idx1=find(idx1==0);
        c_data=B.data(idx1);        
        if flipped
            [c_indices c_partition]=extract_indices(B,[b_outer b_inter]);
            c_indices=c_indices(idx1,:);            
            c_dims=[dims(b_outer) dims(b_inter)];
        else
            [c_indices c_partition]=extract_indices(B,[b_inter b_outer]);
            c_indices=c_indices(idx1,:);            
            c_dims=[dims(b_inter) dims(b_outer)];
        end
    elseif ~isempty(A.outer_idx) %inter/outer product
        
        if flipped
            [c_indices c_partition]=extract_indices(B,[b_outer b_inter b_inter]);            
            c_dims=[dims(b_outer) dims(b_inter) dims(b_inter) ];
        else
            [c_indices c_partition]=extract_indices(B,[b_inter b_inter b_outer ]);     
            c_dims=[dims(b_inter) dims(b_inter) dims(b_outer)];
        end
        c_data=B.data;
    else %inter/inter product
        idx1=B.extract_indices(b_inter(1));
        idx2=B.extract_indices(b_inter(2));
        idx1=idx2-idx1;
        idx1=find(idx1==0);
        c_data=B.data(idx1);
        c_indices=B.indices(idx1,:);
        c_dims=B.dims;
        c_partition=B.partition;
    end
    C=SparseMIA(c_data,c_indices,c_dims,c_partition);
    C=C.consolidate_indices();
elseif ~isempty(A.inner_idx)
    if ~isempty(A.outer_idx) %inner/outer product
        
        if flipped
            [c_indices c_partition]=extract_indices(B,[b_inner b_outer]);
            c_dims=[dims(b_inner) dims(b_outer)];
        else
            [c_indices c_partition]=extract_indices(B,[b_outer b_inner]);
            c_dims=[dims(b_outer) dims(b_inner)];
        end
        C=SparseMIA(B.data,c_indices,c_dims, c_partition);
        C=C.consolidate_indices();
    else %inner/inner product
        
        if flipped
            C=mtimes@SparseMIA(B,A);
        else
            C=mtimes@SparseMIA(A,B);
        end
        
        
    end
else %outer/outer product
    if flipped
        C=mtimes@SparseMIA(B,A);
    else
        C=mtimes@SparseMIA(A,B);
    end
end
