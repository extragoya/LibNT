function error_check_mult(A,B)

a_size=size(A);
if numel(a_size)==2
    if a_size(2)==1;
        a_size=a_size(1);
    end
end
b_size=size(B);
if numel(b_size)==2
    if b_size(2)==1;
        b_size=b_size(1);
    end
end

a_inner_idx=A.inner_idx;
b_inner_idx=B.inner_idx;

if length(a_inner_idx)~=length(b_inner_idx)
    error('Number of indices to undergo an inner product must be the same for each MIA');
end

a_inter_idx=A.inter_idx;
b_inter_idx=B.inter_idx;

if length(a_inter_idx)~=length(b_inter_idx)
    error('Number of indices to undergo an inter product must be the same for each MIA');
end

for i=1:length(a_inner_idx)
    if a_size(a_inner_idx(i))~=b_size(b_inner_idx(i))
        error(['Index ' int2str(a_inner_idx(i)) ' of A must be the same length as index ' int2str(b_inner_idx(i)) ' of B to be inner producted']);
    end
    
end


for i=1:length(a_inter_idx)
    if a_size(a_inter_idx(i))~=b_size(b_inter_idx(i))
        error(['Index ' int2str(a_inter_idx(i)) ' of A must be the same length as index ' int2str(b_inter_idx(i)) ' of B to be inter producted']);
    end
end