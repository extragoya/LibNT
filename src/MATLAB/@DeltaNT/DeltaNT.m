classdef (InferiorClasses = {?DenseNT,?SparseNT}) DeltaNT <DenseNT
    
    
    properties
        
        delta_dim
        
    end
    
    methods
        
        function obj = DeltaNT(varargin)
           if nargin==1 
                arg=varargin{1};
                if isnumeric(arg)
                    
                    obj.delta_dim=arg;
                    
                    
                else
                    error('DeltaNTs can only take one numeric argument within constructor.')
                end
            else
                error('DeltaNTs can only take one numeric argument within constructor.')
            end
            
            
        end
        function ret=size(obj)            
            
            ret=dims(obj);        
            
        end
        function ret=dims(obj)
            ret=[obj.delta_dim obj.delta_dim];       
        end
        
        
        
    
        
        
        
        
        function ret=dimensionality(obj)
            ret=obj.delta_dims^2;
            
        end
        
        
        function B=toLattice(A,row_idx,col_idx,depth_idx)
            
            B=DeltaLattice;
            B.delta_dim=A.delta_dim;
            B.row_idx=numel(row_idx);
            B.col_idx=numel(col_idx);
            B.depth_idx=numel(depth_idx);
        end
       
    end
    
    methods(Static)
        function m_order=order()
            m_order=2;
        end
    
    end
    
   
end