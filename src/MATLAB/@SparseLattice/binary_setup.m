function C = binary_setup(A,B,fhandle,f_intersect)
% binary_setup 
%Original programmed by Dileepan Joseph, updated by Adam Harrison
%UofA Electronic Imaging Lab 2010

%performs checks and if necessary some preprocessing to ensure parameters
%are valid for binary operations and calls either intersect or union

if~exist('f_intersect','var')
    f_intersect=0;
end

[A B IsSparseLattice] = chkbinary(A,B);

if all(IsSparseLattice)
    if isequal(size(A), size(B))
        if f_intersect
            C = intersect(A,B,fhandle);
        else
            C = union(A,B,fhandle);
        end
    else
        error('Operands must be the same size.');
    end
elseif IsSparseLattice(1)
    
    if numel(B) == 1
        C=A;
        C.vals=fhandle(C.vals,B);
        
    elseif isequal(size(A),size(B))
        C = fhandle(full(A),B);
    else
        error('Second operand is not compatible. Valid operands are scalars or matrices/lattices of equal size');
    end
elseif IsSparseLattice(2)
    if numel(A) == 1
        C=B;
        C.vals=fhandle(C.vals,A);
        
    elseif isequal(size(A),size(B))
        C = fhandle(A,full(B));
    else
        
        error('First operand is not compatible. Valid operands are scalars or matrices/lattices of equal size');
    end
else
    error('Neither operand is a sparse lattice.');
end
