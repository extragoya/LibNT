classdef Lattice
    properties (SetAccess = protected)
        vals    %values
        m
        n
        p
        
    end
    
    methods
        
        function obj = Lattice(varargin)
            if nargin==0
                obj.vals=[];
            elseif nargin==1
                arg=varargin{1};
                if isa(arg,'Lattice');
                    obj=arg;
                elseif isnumeric(arg) && numel(size(arg))>1 && numel(size(arg))<4
                    %create a spares lattice based from a full matrix or lattice
                    obj.vals=arg;
                    if numel(size(arg))==2
                        [obj.m obj.n]=size(arg);
                        obj.p=1;
                    else
                        [obj.m obj.n obj.p]=size(arg);
                    end
                else
                    error('Only 2 or 3-indexed numeric arrays may be converted to a Lattice object')
                end
                
            else
                error('Lattice constructor called with incomptabile parameters. Please view help.')
            end
        end % Lattice
        
        
        function ret = subsref(obj,s)
            % SUBSREF Implementing the following syntax:
            % obj(i,j,k)
            dummy=0;
            switch s(1).type
                
                
                case '()'
                    if length(s)==1
                        ret=obj.vals(s(1).subs{:});
                        
                        return
                    else
                        ret = builtin('subsref',obj,s);
                    end
                case '.'
                    
                    ret = builtin('subsref',obj,s);
                    
                otherwise
                    error('Specify value for x as obj(x)')
            end
        end % subsref
        
        function ret=size(obj)
            
            ret=size(obj.vals);
            
        end
        
    end
    methods
        C=plus(A,B)
        C=minus(A,B)
        C=mtimes(A,B)
        B=mldivide(A,C);
        A=transpose(A);
        A=ctranspose(A);
    end
    
end %classdef