function A=consolidate_indices(A)

i=1;

dims=A.dims;

new_partition=zeros(1,length(dims));
lin_length=1;
k=1;
for i=1:length(dims)
    new_length=lin_length*A.dims(i);
    if new_length>intmax('int64')
        lin_length=A.dims(i);
        new_partition(k)=i-1;
        k=k+1;
    else
        lin_length=new_length;
    end
end
new_partition(k)=i;
new_partition=new_partition(1:k);

new_indices=int64(zeros(size(A.indices,1),numel(new_partition)));

start_index=1;
for i=1:numel(new_partition)
    running_product=1;
    for j=start_index:new_partition(i)
        if j==start_index
            new_indices(:,i)=A.pull_index(j);
        else
            running_product=running_product*dims(j-1);
            new_indices(:,i)=new_indices(:,i)+running_product.*(A.pull_index(j)-1);
        end
    
    end
    start_index=new_partition(i)+1;
    
end
A.indices=new_indices;
A.partition=new_partition;