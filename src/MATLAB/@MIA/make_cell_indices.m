function split_indices=make_cell_indices(obj,indices)
%splits a char array, e.g., 'i!kj' into a cell array, e.g., {i},{!k},{j}
%and also throws an error with invalid input

if(ischar(indices))
    %check that the number of indices is the same as the MIA order
    
    reg='!*[a-z]';
    not_reg='(?=[^a-z])[^!]';
    split_indices=regexp(indices,reg,'match');
    invalid_indices=regexp(indices,not_reg,'match');
    if ~isempty(invalid_indices)
        error('Only input lower-case letters and ''!''') 
    end
    if length(split_indices)~=obj.order
       error('Number of indices must be the same as MIA order') 
    end
    
    
else
    error('Must index MIA with a char array, e.g., ''ijk''')
end