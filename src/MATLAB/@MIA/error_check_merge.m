function error_check_merge(A,B)

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

a_merge_idx=A.merge_idx;
b_merge_idx=B.merge_idx;

if length(a_merge_idx)~=length(b_merge_idx)
    error('Number of dimensions must be the same for each MIA being merged');
end


for i=1:length(b_merge_idx)
    if a_size(a_merge_idx(i))~=b_size(b_merge_idx(i))
        error(['For merger, index ' int2str(a_merge_idx(i)) ' of A must be the same length as index ' int2str(b_merge_idx(i))]);
    end
end