function C=plus(A,B)

permute_idx=pull_merge_indices(A,B);
error_check_merge(A,B,permute_idx);
C_mia=A.m_mia.do_plus(B.m_mia,permute_idx);
C=NTExpr(C_mia,A.m_indices);

