function ret = delete(A,varargin)
%given list of indices, selects the sub NT corresponding to the list. All
%selected indices are collapsed to one index. E.g., if a fourth-degree NT is
%selected with A.select(xs,ys,:,:), this will return a third-degree NT,
%where the first two indices, are collapsed into a length(xs) index,
%corresponding to A(xs(1),ys(1),:,:) to A(xs(n),ys(n),:,:);
ret_compl=0;
[selecting_indices, indices_to_select]=check_selection_args(A,varargin{:});
no_number=length(selecting_indices);
all_indices=1:A.order;
other_indices=setdiff(all_indices,selecting_indices); %get the indices not being selected
A=A.permute([other_indices selecting_indices]); %permute A so that the selecting indices are in highest significance
a_dims=A.dims;

selecting_dims=zeros(1,no_number);
selecting_dims(1)=1;
for i=1:no_number-1
    selecting_dims(i+1)=a_dims(selecting_indices(i));
end
divisor=1; %divisor to use to isolate selecting indices, and remaining indices
other_dims=zeros(1,length(other_indices));
for i=1:length(other_indices)
    divisor=divisor*a_dims(other_indices(i));
    other_dims(i)=a_dims(other_indices(i)); %computer other dims
end
divisor=int64(divisor);

indices_to_select=indices_to_select-1;
indices_to_select=int64(indices_to_select*selecting_dims'); %get the linearized version of the indices
new_inds=A.indices-1;
temp_inds=idivide(new_inds,divisor); %isolate only the selecting indices

[LIA] = ismember(temp_inds,indices_to_select); %see which indices in A match with the selecting indices
new_inds=new_inds+1;
new_inds=new_inds(~LIA); %pull non deleted indices
new_data=A.data(~LIA);   %pull non deleted data

ret=SparseNT(new_data,new_inds,A.dims(),[other_indices selecting_indices],true);


end

