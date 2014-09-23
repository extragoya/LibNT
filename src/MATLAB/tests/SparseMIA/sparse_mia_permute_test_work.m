function [A, A2]=sparse_mia_permute_test_work(newLinIdx,dims,A,A2)


if nargin==2
    a_data=rand(dims);
    a_data(a_data<0.5)=0;
    
    DenseA=MIA(a_data);
    A=SparseMIA(DenseA);
    A2=SparseMIA(DenseA);
    
end
A=A.permute(newLinIdx);
%this test assumes that changeLinIdx and sort are working correctly
A2=A2.changeLexOrder(newLinIdx);
A2=A2.sort();

   











