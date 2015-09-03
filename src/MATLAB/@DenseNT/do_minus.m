function C= do_minus(A,B,permute_idx)
%perform addition based on the permutation idx of B
if(~isa(B,'DenseNT'))
    error('Unsupported class in do_minus');
end
if(isa(B,'SparseNT')) %if B is a sparseNT we can safely convert it to sparseNT because it's the same size as A
    B=DenseNT(B);
end

do_permute=true;
a_size=size(A);
if numel(a_size)==2
    if a_size(2)==1;
        do_permute=false;
    end
end

if do_permute
    
    B_data=permute(B.data,permute_idx);
end


C_data=A.data-B_data;
C=DenseNT(C_data);
end

