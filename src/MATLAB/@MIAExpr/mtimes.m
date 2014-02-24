function CExpr = mtimes(A,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[a_inner_idx, a_inter_idx, a_outer_idx, b_inner_idx, b_inter_idx, b_outer_idx]=pull_mult_indices(A,B);
error_check_mult(A,B,a_inner_idx, a_inter_idx, b_inner_idx, b_inter_idx) %make sure indices match up properly for multiplication



A_mia=A.m_mia;
B_mia=B.m_mia;
a_size=size(A_mia);
b_size=size(B_mia);



%extract index ranges for A and B. Also calculate ranges after A and B are
%flattened into lattices
a_inter_size=a_size(a_inter_idx);
a_outer_size=a_size(a_outer_idx);
b_outer_size=b_size(b_outer_idx);


A_lat=A_mia.toLattice(a_outer_idx, a_inner_idx, a_inter_idx);
B_lat=B_mia.toLattice(b_inner_idx, b_outer_idx, b_inter_idx);



%perform lattice product
C_lat=A_lat*B_lat;
clear A_lat B_lat
%convert result back to an MIA
c_indices=A.m_indices(a_outer_idx);
c_indices=[c_indices B.m_indices(b_outer_idx)];
c_indices=[c_indices A.m_indices(a_inter_idx)];

C=C_lat.toMIA(a_outer_size, b_outer_size, a_inter_size);
%create a MIAExpr using resulting MIA and indices
CExpr=MIAExpr(C,c_indices);


end

