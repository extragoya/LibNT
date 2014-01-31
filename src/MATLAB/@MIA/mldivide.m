function B=mldivide(A,C)
%MIA inversion 
%
%B=A\C computes the solution of A*B=C. Indices of A undergoing inner, inter, and
%outer product must be set by member function set_idx. If indices are not 
%set, A and B undergo by default an inter product for each corresponding 
%element.
%
%Order of indices after multiplication are [A.outer_idx B.outer_idx
%inter_idx]
%
%Developed by Adam Harrison
%U of A Electronic Imaging Lab, 2011


error_check_mldivide(A,C); %make sure indices match up properly for multiplication

a_inner_idx=A.inner_idx;
c_inner_idx=C.inner_idx;

a_size=size(A);
c_size=size(C);

a_inter_idx=A.inter_idx;
c_inter_idx=C.inter_idx;

a_outer_idx=A.outer_idx;
c_outer_idx=C.outer_idx;

%extract index ranges for A and B. Also calculate ranges after A and B are
%flattened into lattices
a_inner_size=a_size(a_inner_idx);
a_inner_length=prod(a_inner_size);
c_inner_size=c_size(c_inner_idx);
c_inner_length=prod(c_inner_size);

a_inter_size=a_size(a_inter_idx);
a_inter_length=prod(a_inter_size);
c_inter_size=c_size(c_inter_idx);
c_inter_length=prod(c_inter_size);

a_outer_size=a_size(a_outer_idx);
a_outer_length=prod(a_outer_size);
c_outer_size=c_size(c_outer_idx);
c_outer_length=prod(c_outer_size);

A_lat=A.toLattice(a_outer_idx, a_inner_idx, a_inter_idx,[a_outer_length,a_inner_length,a_inter_length]);
C_lat=C.toLattice(c_inner_idx, c_outer_idx, c_inter_idx,[c_inner_length c_outer_length c_inter_length]);

clear A C

%perform lattice inversion
B_lat=A_lat\C_lat;
clear A_lat C_lat

%convert result back to an MIA
B_data=B_lat.vals;
B_data=squeeze(B_data);
new_dims=[a_inner_size c_outer_size a_inter_size];
if length(new_dims)>1
    B_data=reshape(B_data,new_dims);
end
B_data=squeeze(B_data);
B=MIA(B_data);

end