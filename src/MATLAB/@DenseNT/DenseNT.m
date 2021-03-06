classdef DenseNT
    
    
    properties
        
        data
        solveInfo=0; %no Info
        
    end
    
    methods
        
        function obj = DenseNT(varargin)
            if nargin==0
                obj.data=[];
                
            elseif nargin==1 || nargin==2
                arg=varargin{1};
                if isnumeric(arg)
                    
                    if nargin==1
                        obj.solveInfo=0;
                    else
                        obj.solveInfo=varargin{2};
                    end
                    
                    obj.data=full(arg);
                    obj.data=squeeze(obj.data);
                    s=size(obj.data);
                    l=length(s);
                    if l==2
                        idx=find(s==1);
                        if idx
                            if idx==1
                                obj.data=obj.data';
                            end
                       
                        end
                    end
                    
                elseif isa(arg,'SparseNT');
                    arg=arg.permute([1:arg.order]);
                    dims=size(arg);
                    obj.data=zeros(dims);
                    indices=arg.indices;      
                    if size(indices,2)>1
                        error('SparseNT dimensions are too large to allow dense representation')
                    end
                                  
                    obj.data(indices)=arg.data;
                    obj.data=squeeze(obj.data);
                    
                    
                    
                else
                    error('Input a numeric array into DenseNT constructor')
                end
            else
                error('DenseNT constructor called with incomptabile parameters. Please view help.')
            end
            
            
        end
        function ret=size(obj)            
            
            ret=size(obj.data);            
            
        end
        function ret=dims(obj)
            ret=size(obj.data);
            if (ret(end)==1)
                ret=ret(1:end-1);
            end
        end
        
        function obj=permute(obj,order)
            obj.data=permute(obj.data,order);
        
        end
        
        function sref = subsref(obj,s)
            switch s(1).type
                case '{}'
                    sref=obj.data(s(1).subs{:});
                case '()'
                    sref=obj.make_expr(s(1).subs{:});
                case '.'
                    sref = builtin('subsref',obj,s);
                otherwise
                    error('Invalid indexing');
                        
            end
        end
        
        function obj = subsasgn(obj,s,val)
            switch s(1).type
                case '()'
                    obj=assign_expr(obj,s(1).subs{:},val);                
                case '.'
                    obj = builtin('subsasgn',obj,s,val);
                otherwise
                    error('Invalid indexing');
                        
            end
        end
        
    
        
        function ret=uminus(obj)
            ret=obj;
            ret.data=-(ret.data);
        end
        
        function m_order=order(obj)
            m_order=numel(size(obj.data));
            %1D arrays will still have a second dimension of 1, so check
            %that
            if m_order==2
                if size(obj.data,2)==1
                    m_order=m_order-1;
                end
            end
        end
        
        function ret=dimensionality(obj)
            ret=prod(obj.dims);
            
        end
        
        isequal=eq(a,b);
        B=toLattice(A,row_idx,col_idx,depth_idx);
        C=do_plus(A,B,permute_idx);
        C=do_minus(A,B,permute_idx);
        C=mtimes(A,B); %just scalar multiplication
        C=mrdivide(A,B); %just scalar division
        B=sqrt(A);
        B=exp(A);
        Amatrix=flatten(A,row_inds,col_inds);
    end
    
    methods (Access=protected)
        Expr=make_expr(A,indices);
        obj=assign_expr(obj,indices,Expr);
        is_equal=compare(A,B,func);
        obj=assign(obj,otherMIA,assign_order);
        
    end
    methods(Access=protected, Static)
        split_indices=make_cell_indices(indices);
        
    end
    
end