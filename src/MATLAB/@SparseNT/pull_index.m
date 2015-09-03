function indices=pull_index(A,idx)
%returns expanded indices of non-zeros indices of A
indices=int64(zeros(size(A.indices,1),length(idx)));
dim=A.dims();
dims=dim(A.lexOrder);
for i=1:length(idx)
    j=1;  
        
   
    temp=int64(A.indices(:,j)-1);
    if (idx(i)>1)
        running_product=prod(dims(1:idx(i)-1));
        
        temp=idivide(temp,int64(running_product));
    end
    temp=mod(temp,dims(idx(i)))+1;
    
    
    indices(:,i)=temp;
    
        
end