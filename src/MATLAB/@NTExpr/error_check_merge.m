function error_check_merge(A,B, permute_idx)


a_size=size(A.m_mia);
if numel(a_size)==2
    if a_size(2)==1;
        a_size=a_size(1);
    end
end
b_size=size(B.m_mia);
if numel(b_size)==2
    if b_size(2)==1;
        b_size=b_size(1);
    end
end




for i=1:length(permute_idx)
    if a_size(i)~=b_size(permute_idx(i))
        error(['For merger, index ' int2str(a_size(i)) ' of A must be the same length as index ' int2str(permute_idx(i))]);
    end
end