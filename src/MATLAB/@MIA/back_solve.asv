function c=back_solve(A)



dims=size(A);
dim=numel(dims)-1;
elements=dims(1)-1;
c=zeros(elements,dim);

%get roots of last row polynomial

idx=zeros(dim+1,dim+1);
baseline=ones(1,dim+1)*elements;
for j=1:dim
    idx(j,:)=baseline;
    baseline(dim+1-j+1)=elements+1;
end
idx(dim+1,:)=ones(1,dim+1)*(elements+1);
idx(dim+1,1)=elements;
idx_cell=cell(1,dim+1);
for j=1:dim+1
    idx_cell{j}=idx(j,:);
end
idx=sub2ind(dims,idx_cell{:});
coeffs=A(idx);
r = roots(coeffs);
c(elements,:)=r;    

contractor=zeros(elements-1,elements);
contractor(1:elements-1,1:elements-1)=eye(elements-1,elements-1);
CR=MIA(contractor);
CR=CR.set_idx(2,[],1);
A.set_idx(1,[],2:dim+1);
A=CR*A;

contractor=zeros(elements,elements+1);
contractor(1:elements,1:elements)=eye(elements,elements);
contractor(elements,elements+1)=1;
contractor=repmat(contractor,[elements elements+1 length(r)]);
contractor=permute(contractor,[3 1 2]);
for i=1:length(r)
    contractor(i,elements,elements)=r(i);
    
end

CR=MIA(contractor);
CR=CR.set_idx(3,1,2);

A=A.cascade_solution(CR);

for i=elements-1:1
    idx=zeros(dim+1,dim+1);
    baseline=ones(1,dim+1)*element;
    for j=1:dim
        idx(j,:)=baseline;
        baseline(dim+1-j+1)=element+1;
    end
    idx(dim+1,:)=ones(1,dim+1)*element+1;    
    idx(dim+1,1)=element;
    idx_cell=cell(1,dim+1);
    for j=1:dim+1
        idx_cell{j}=idx(j,:);
    end
    idx=sub2ind(dims,idx_cell{:});
    coeffs=A(idx);
    r = roots(coeffs);
    
end