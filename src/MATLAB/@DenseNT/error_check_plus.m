function error_check_plus(A,B)

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

a_plus_idx=A.a_plus_idx;
b_plus_idx=B.b_plus_idx;

if length(a_plus_idx)~=length(b_plus_idx)
    error('Number of dimensions must be the same for each NT being added');
end


for i=1:length(a_plus_idx)
    if a_size(a_plus_idx(i))~=b_size(b_plus_idx(i))
        error(['For addition, index ' int2str(a_plus_idx(i)) ' of A must be the same length as index ' int2str(b_plus_idx(i))]);
    end
end