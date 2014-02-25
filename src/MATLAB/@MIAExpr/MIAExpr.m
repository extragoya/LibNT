classdef MIAExpr
    
    properties
    
        m_mia %references MIA that is part of expression
        m_indices %the indices used in the expression
    
    end
    
    methods
        
         function obj = MIAExpr(varargin)
            if nargin==2
                arg1=varargin{1};
                arg2=varargin{2};
                if iscell(arg2) && isa(arg1,'MIA')                    
                    obj.m_mia=arg1;
                    obj.m_indices=arg2;
                else
                    error('An MIA and cell array of strings are the only valid arguments to the MIAExpr class')
                end
                
            else
                error('Only one MIA and string argument can be inputted as index expressions, e.g., ''ijk''')
            end
         end
         
         C=mtimes(A,B);
         C=plus(A,B);
         C=minus(A,B);
    end
    methods (Access=protected)
        [a_inner_idx, a_inter_idx, a_outer_idx, b_inner_idx, b_inter_idx, b_outer_idx]=pull_mult_indices(A,B);
        error_check_mult(A,B,a_inner_idx, a_inter_idx, b_inner_idx, b_inter_idx)
        error_check_mldivide(A,C);
        error_check_merge(A,B,permute_idx);
        error_check_flatten(A,row_idx,col_idx);
    end
end