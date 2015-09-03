function C=minus(A,B)

permute_idx=pull_merge_indices(A,B);
error_check_merge(A,B,permute_idx);
C_NT=A.m_mia.do_minus(B.m_mia,permute_idx);
C=NTExpr(C_NT,A.m_indices);

