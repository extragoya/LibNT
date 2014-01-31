function C=plus(A,B)


error_check_merge(A,B);

B_data=B.data;
A_data=A.data;

do_permute=true;
a_size=size(A);
if numel(a_size)==2
    if a_size(2)==1;
        do_permute=false;
    end
end

if do_permute
    A_data=permute(A_data,A.merge_idx);
    B_data=permute(B_data,B.merge_idx);
end

C_data=A_data+B_data;
C=MIA(C_data);