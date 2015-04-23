function [B, DenseB]=sparse_mia_assign_test_work(a_dims,a_indices,b_indices,lexOrder)
%a_indices, b_indices, the string of indices, e.g. 'ijk'. permute_idx, the
%actual numerical indices by which A and B are matched. That is
%A{1,2,3}=B{permute_idx(1),permute_idx(2),permute_idx(3)};

if(nargin==3)
    lexOrder=1:length(a_dims);
end
a_data=rand(a_dims);
a_data(a_data<0.5)=0;
DenseA=MIA(a_data);
A=SparseMIA(DenseA);
A=A.permute(lexOrder);
B=SparseMIA;
DenseB=MIA;
B(b_indices)=A(a_indices);
DenseB(b_indices)=DenseA(a_indices);





