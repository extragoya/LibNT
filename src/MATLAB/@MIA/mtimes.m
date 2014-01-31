function C=mtimes(A,B)
%MIA multiplication 
%
%C=A*B multiplies two MIAs together. Indices undergoing inner, inter, and
%outer product must be set by member function set_idx. If indices are not 
%set, A and B undergo by default an inter product for each corresponding 
%element.
%
%Order of indices after multiplication are [A.outer_idx inter_idx B.outer_idx]
%
%Developed by Adam Harrison
%U of A Electronic Imaging Lab, 2011

if ~isa(A,'MIA') || ~isa(B,'MIA')
    error('MIAs can only be multipled with classes or subclasses of MIA');
end

error_check_mult(A,B); %make sure indices match up properly for multiplication

a_inner_idx=A.inner_idx;
b_inner_idx=B.inner_idx;

a_size=size(A);
b_size=size(B);

a_inter_idx=A.inter_idx;
b_inter_idx=B.inter_idx;

a_outer_idx=A.outer_idx;
b_outer_idx=B.outer_idx;

%extract index ranges for A and B. Also calculate ranges after A and B are
%flattened into lattices
a_inner_size=a_size(a_inner_idx);
a_inner_length=prod(a_inner_size);
b_inner_size=b_size(b_inner_idx);
b_inner_length=prod(b_inner_size);

a_inter_size=a_size(a_inter_idx);
a_inter_length=prod(a_inter_size);
b_inter_size=b_size(b_inter_idx);
b_inter_length=prod(b_inter_size);

a_outer_size=a_size(a_outer_idx);
a_outer_length=prod(a_outer_size);
b_outer_size=b_size(b_outer_idx);
b_outer_length=prod(b_outer_size);

A_lat=A.toLattice(a_outer_idx, a_inner_idx, a_inter_idx,[a_outer_length,a_inner_length,a_inter_length]);
B_lat=B.toLattice(b_inner_idx, b_outer_idx, b_inter_idx,[b_inner_length b_outer_length b_inter_length]);

clear A B

%perform lattice product
C_lat=A_lat*B_lat;
clear A_lat B_lat

%convert result back to an MIA
C_data=C_lat.vals;
C_data=permute(C_data,[1 3 2]);
C_data=squeeze(C_data);
new_dims=[a_outer_size a_inter_size b_outer_size];
if length(new_dims)>1
    C_data=reshape(C_data,new_dims);
end
C_data=squeeze(C_data);
C=MIA(C_data);

end