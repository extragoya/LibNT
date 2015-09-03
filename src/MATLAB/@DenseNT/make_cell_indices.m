function split_indices=make_cell_indices(indices)
%splits a char array, e.g., 'i!kj' into a cell array, e.g., {i},{!k},{j}
%and also throws an error with invalid input

if(ischar(indices))
    %check that the number of indices is the same as the NT order
    
    %reg='(!~,)*[a-z]';
    
    not_reg='(?=[^a-z])[^!,]';
    split_indices=regexp(indices,',','split');
    invalid_indices=regexp(indices,not_reg,'match');
    if ~isempty(invalid_indices)
        error('Only input lower-case letters and ''!'' or '',''') 
    end
    letter_start=regexp(split_indices,'[a-z]', 'once');
    test=cellfun('isempty',letter_start);
    if sum(test)>0
        error('Must input at least one letter as an NT index');
    end
    
    letter_start=regexp(split_indices,'[a-z]!', 'once');
    test=cellfun('isempty',letter_start);
    if sum(test)~=length(split_indices)
        error('Any exclamation points must come before letters in an NT index');
    end
    
    
    
    
else
    error('Must index NT with a char array, e.g., ''ijk''')
end