function [ret, ret_compl ] = select_subNT(A,varargin)
%given list of indices, selects the sub NT corresponding to the list. All
%selected indices are collapsed to one index. E.g., if a fourth-degree NT is
%selected with A.select(xs,ys,:,:), this will return a third-degree NT,
%where the first two indices, are collapsed into a length(xs) index,
%corresponding to A(xs(1),ys(1),:,:) to A(xs(n),ys(n),:,:);
ret_compl=0;
[selecting_indices, indices_to_select]=check_selection_args(A,varargin);
no_number=length(selecting_indices);
all_indices=1:A.order;
other_indices=setdiff(all_indices,selecting_indices); %get the indices not being selected
A=A.permute([other_indices selecting_indices]); %permute A so that the selecting indices are in highest significance
a_dims=A.dims;
no_number=no_number-1;
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

[LIA, LIB] = ismember(temp_inds,indices_to_select); %see which indices in A match with the selecting indices
LIB=int64(LIB);
new_inds=new_inds(LIA); %pull only the matching indices
new_inds=mod(new_inds,divisor); %pulled indices should only correspond to remaining indinces, ie, not the selecting indices
new_inds=new_inds+(LIB(LIA)-1)*divisor; %add the new selecting indices dimension
new_inds=new_inds+1;
new_data=A.data(LIA);   %pull only corresponding data
new_dims=[length(indices_to_select) other_dims]; %new set of dimensions
new_order=length(new_dims);
ret=SparseNT(new_data,new_inds,new_dims,[2:new_order 1],true);

if nargout==2
    selecting_dims=prod(selecting_dims)*a_dims(selecting_indices(end));
    if(selecting_dims<=length(temp_inds)) %only perform this if selecting dims is hypersparse
        compl_indices=true(selecting_dims,1);
        compl_indices(indices_to_select+1)=false;
        compl_indices_to_select=1:selecting_dims;
        compl_indices_to_select(compl_indices)=1:selecting_dims-length(indices_to_select);%compl_indices_to_select(compl_indices);
        compl_indices_to_select=compl_indices_to_select-1;
        temp_inds=temp_inds(~LIA);
        temp_inds=int64(compl_indices_to_select(temp_inds));       
        temp_inds=temp_inds';
    else
        temp_inds=temp_inds(~LIA);
        for i=length(indices_to_select):-1:1
            temp_idx=temp_inds>indices_to_select(i);
            temp_inds(temp_idx)=temp_inds(temp_idx)-1;
        end
        
    end
    new_inds=A.indices-1;
    new_inds=new_inds(~LIA); %pull only the matching indices
    new_inds=mod(new_inds,divisor); %pulled indices should only correspond to remaining indinces, ie, not the selecting indices
    new_inds=new_inds+temp_inds*divisor; %add the new selecting indices dimension
    new_inds=new_inds+1;
    new_data=A.data(~LIA);   %pull only corresponding data
    new_dims=[selecting_dims-length(indices_to_select) other_dims]; %new set of dimensions    
    ret_compl=SparseNT(new_data,new_inds,new_dims,[2:new_order 1],true);
end

end

