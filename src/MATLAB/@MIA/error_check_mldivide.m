function error_check_mldivide(A,C)

a_size=size(A);
if numel(a_size)==2
    if a_size(2)==1;
        a_size=a_size(1);
    end
end
c_size=size(C);
if numel(c_size)==2
    if c_size(2)==1;
        c_size=c_size(1);
    end
end


a_inter_idx=A.inter_idx;
c_inter_idx=C.inter_idx;

if length(a_inter_idx)~=length(c_inter_idx)
    error('Number of indices to undergo an inter product must be the same for each MIA');
end


for i=1:length(a_inter_idx)
    if a_size(a_inter_idx(i))~=c_size(c_inter_idx(i))
        error(['Index ' int2str(a_inter_idx(i)) ' of A must be the same length as index ' int2str(c_inter_idx(i)) ' of B to be inter producted']);
    end
end