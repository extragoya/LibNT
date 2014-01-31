function [A,B,IsSparseLattice] = chkbinary(A,B)
% CHKBINARY (check operands for sparse lattices and convert is possible)
% Originally programmed by Dileepan Joseph
% Updated by Adam Harrison
% UofA Electronic Imaging Lab 2010

IsSparseLattice = logical([0 0]);

switch class(A)
case 'SparseLattice'
   IsSparseLattice(1) = 1;
case 'sparse'
   A = sparray(A);
   IsSparseLattice(1) = 1;
otherwise
   IsSparseLattice(1) = 0;
end

switch class(B)
case 'SparseLattice'
   IsSparseLattice(2) = 1;
case 'sparse'
   B = sparray(B);
   IsSparseLattice(2) = 1;
otherwise
   IsSparseLattice(2) = 0;
end
