function [A, A2]=sparse_mia_changelinIdx_test_work(newLinIdx,dims,A,A2)

if nargin==2
    a_data=rand(dims);
    a_data(a_data<0.5)=0;
    
    DenseA=MIA(a_data);
    A=SparseMIA(DenseA);
    A2=SparseMIA(DenseA);
    
end
A=A.changeLexOrder(newLinIdx);
new_dims=dims(newLinIdx);

inds=A2.pull_index(1:numel(dims));
inds(:,A2.lexOrder)=inds; %get indices in defaul lexOrder
inds=inds(:,newLinIdx);
inds=inds';
inds = inds-1;
A2.indices = (int64(cumprod([1 new_dims(1:end-1)])*double(inds))+1)';
A2.lexOrder=newLinIdx;

   











