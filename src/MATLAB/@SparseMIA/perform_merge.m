function C=perform_merge(A,B,fhandle)

error_check_merge(A,B);
[A.indices A.partition]=A.extract_indices(A.merge_idx);
[A.indices idx]=sortrows(A.indices);
A.data=A.data(idx);
A.dims=A.dims(A.merge_idx);
[B.indices B.partition]=B.extract_indices(B.merge_idx);
[B.indices idx]=sortrows(B.indices);
B.data=B.data(idx);
B.dims=B.dims(B.merge_idx);
C=union(A,B,fhandle);


