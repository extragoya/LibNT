function indices=pull_index(A,idx)

indices=int64(zeros(size(A.indices,1),length(idx)));

for i=1:length(idx)
    j=1;
    while idx(i)>A.partition(j)
        j=j+1;
        
    end
    if j>1
        start=A.partition(j-1)+1;
    else
        start=1;
    end
    temp=int64(A.indices(:,j)-1);
    if (idx(i)>1)
        running_product=prod(A.dims(start:idx(i)-1));
        
        temp=idivide(temp,int64(running_product));
    end
    temp=mod(temp,A.dims(idx(i)))+1;
    
    
    indices(:,i)=temp;
    
        
end