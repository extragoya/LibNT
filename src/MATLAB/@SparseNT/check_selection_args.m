function [selecting_indices, indices_to_select]=check_selection_args(A,varargin)
%given list of indices, selects the sub NT corresponding to the list. All
%selected indices are collapsed to one index. E.g., if a fourth-degree NT is
%selected with A.select(xs,ys,:,:), this will return a third-degree NT,
%where the first two indices, are collapsed into a length(xs) index,
%corresponding to A(xs(1),ys(1),:,:) to A(xs(n),ys(n),:,:);

if length(varargin)~=A.order
    error('Your selection must have as many entries as the NT order');
end
no_number=0;
selecting_indices=[];
indices_to_select=[];
for i=1:A.order
    if isnumeric(varargin{i})
        no_number=no_number+1;
        if(numel(size(varargin{i}))~=2)
            error('Only input vector as indices');
        end
        if (size(varargin{i},2)~=1)
            varargin{i}=varargin{i}';
        end
        if (size(varargin{i},2)~=1 )
            error('Only input vector as indices');
        end
        if (isscalar(varargin{i}) && varargin{i}==-1)
            continue;
        end
        length_indices=size(varargin{i},1);
        
    elseif ischar(varargin{i}) && (varargin{i}=='a' || varargin{i}==':')
    else
        error('Only input numeric indices or '':'' when selecting an NT.')
    end
        
end

if ~no_number %if all ':', we can just return and do nothing
    return;
end
indices_to_select=zeros(length_indices,no_number);
selecting_indices=zeros(1,no_number);

no_number=1;
for i=1:A.order
    if ~isnumeric(varargin{i}) 
        continue;
    end
    if size(varargin{i},1)~=length_indices
        error('Input vectors of same size to select indices');
    end
    selecting_indices(no_number)=i;
    indices_to_select(:,no_number)=varargin{i};
    no_number=no_number+1;
end


end

