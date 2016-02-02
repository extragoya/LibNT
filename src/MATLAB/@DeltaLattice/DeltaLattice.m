classdef (InferiorClasses = {?Lattice,?SparseLattice}) DeltaLattice
    
    
    properties
        
        delta_dim
        row_idx
        col_idx 
        depth_idx
    end
    
    methods
        
        function obj = DeltaLattice(varargin)
           if nargin==1 
                arg=varargin{1};
                if isnumeric(arg)
                    
                    obj.delta_dim=arg;
                    
                    
                else
                    error('DeltaLattice can only take one numeric argument within constructor.')
                end
           elseif nargin==0
               obj.delta_dim=0;
           else
                error('DeltaLattice can only take one numeric argument within constructor.')
            end
            
            
        end
        
        
        C=mtimes(A,B);
        function ret=size(obj)            
            m=obj.delta_dim^(obj.row_idx);
            n=obj.delta_dim^(obj.col_idx);
            p=obj.delta_dim^(obj.depth_idx);
            ret=[m n p];        
            
        end
        function ret=p(obj)
            ret= obj.delta_dim;
        end
    end
    
    methods(Access=protected)
        C=sparse_times(A,B);
    end
end