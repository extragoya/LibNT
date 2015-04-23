function error_check_mldivide(A,B,a_inter_idx,b_inter_idx)

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




if length(a_inter_idx)~=length(b_inter_idx)
    error('Number of indices to undergo an inter product must be the same for each MIA');
end


for i=1:length(a_inter_idx)
    if a_size(a_inter_idx(i))~=b_size(b_inter_idx(i))
        error(['Index ' int2str(a_inter_idx(i)) ' of A must be the same length as index ' int2str(b_inter_idx(i)) ' of B for this solution of equations']);
    end
end