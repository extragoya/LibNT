function [indices new_partition]=extract_indices(A,idx,linear)

if ~exist('linear','var')
    linear=true;
end

if linear
    new_partition=zeros(1,length(idx));
    lin_length=1;
    k=1;
    for i=1:length(idx)
        new_length=lin_length*A.dims(idx(i));
        if new_length>intmax('int64')
            lin_length=A.dims(idx(i));
            new_partition(k)=i-1;
            k=k+1;
        else
            lin_length=new_length;
        end
    end
    new_partition(k)=i;
    new_partition=new_partition(1:k);
    
    indices=int64(zeros(size(A.indices,1),k));
    
    start_idx=1;
    for i=1:k
        running_product=1;
        for j=start_idx:new_partition(k)
            t_indice=A.pull_index(idx(j));
            if j==start_idx
                indices(:,k)=t_indice;
            else
                t_indice=t_indice-1;
                indices(:,k)=indices(:,k)+t_indice*running_product;
            end
            running_product=running_product*A.dims(idx(j));
        end
        start_idx=new_partition(k)+1;
    end
else
    new_partition=[];
    indices=int64(zeros(size(A.indices,1),length(idx)));
    for i=1:length(idx)
        indices(:,i)=A.pull_index(idx(i));
    end
end
