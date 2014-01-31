classdef MIA
    
    
    properties
        
        data
        inner_idx
        inter_idx
        outer_idx
        merge_idx
    end
    
    methods
        
        function obj = MIA(varargin)
            if nargin==0
                obj.data=[];
                obj.inner_idx=[];
                obj.outer_idx=[];
                obj.inter_idx=[];
            elseif nargin==1
                arg=varargin{1};
                if isnumeric(arg)
                    obj.data=arg;
                    obj.data=squeeze(obj.data);
                    s=size(obj.data);
                    l=length(s);
                    if l==2
                        idx=find(s==1);
                        if idx
                            if idx==1
                                obj.data=obj.data';
                            end
                        l=l-1;
                        end
                    end
                    obj.merge_idx=1:l;
                    obj.inter_idx=1:l;
                    obj.inner_idx=[];
                    obj.outer_idx=[];
                elseif isa(arg,'SparseMIA');
                    dims=arg.dims;
                    obj.data=zeros(dims);
                    indices=arg.indices;      
                    if size(indices,2)>1
                        error('SparseMIA dimensions are too large to allow dense representation')
                    end
                                  
                    obj.data(indices)=arg.data;
                    obj.data=squeeze(obj.data);
                    s=size(obj.data);
                    l=length(s);
                    obj.inter_idx=1:l;
                    obj.inner_idx=[];
                    obj.outer_idx=[];
                else
                    error('Input a numeric array into MIA constructor')
                end
            else
                error('MIA constructor called with incomptabile parameters. Please view help.')
            end
            
            
        end
        function ret=size(obj)            
            
            ret=size(obj.data);            
            
        end
        function obj=set_idx(obj,inner,inter,outer)
            
            if ~exist('inner','var')
                inner=[];
            end
            if ~exist('inter','var')
                inter=[];
            end
            if ~exist('outer','var')
                outer=[];
            end
            range=size(obj);
            dim=numel(range);
            if dim==2
                if find(range==1)
                    dim=1;
                end
                
            end
            if (~isnumeric(inner) || ~isvector(inner))&&~isempty(inner)
                error('Input numeric vector arrays or empty brackets to specify inner product indices');               
            end
            if (~isnumeric(inter) || ~isvector(inter))&&~isempty(inter)
                error('Input numeric vector arrays or empty brackets to specify inter product indices');               
            end
            if (~isnumeric(outer) || ~isvector(outer))&&~isempty(outer)
                error('Input numeric vector arrays or empty brackets to specify outer product indices');               
            end
            
            if sum(floor(inner)-inner)~=0 || sum(floor(inter)-inter)~=0 || ...
                    sum(floor(outer)-outer)~=0
                error('Indices must be integer valued')
            end
            if ~isempty(inner)&&(min(inner)<1 || max(inner)>dim)
                error(['Inner product indices must be between 1 and ' int2str(dim)])
            end
            if ~isempty(inter)&&(min(inter)<1 || max(inter)>dim)
                error(['Inner product indices must be between 1 and ' int2str(dim)])
            end
            if ~isempty(outer)&&(min(outer)<1 || max(outer)>dim)
                error(['Inner product indices must be between 1 and ' int2str(dim)])
            end
            
            if ~isempty(intersect(inner,inter))
                error('Inner and inter product have repeated index')
            end
            if ~isempty(intersect(inner,outer))
                error('Inner and outer product have repeated index')
            end
            if ~isempty(intersect(inter,outer))
                error('Inter and outer product have repeated index')
            end
            indices=1:dim;
            indices(inner)=0;            
            indices(inter)=0;
            indices=find(indices);
            if ~isempty(setdiff(indices, outer)) && ~isempty(outer)
                error('Not all indices have been set');
            end
            obj.inner_idx=inner;
            obj.inter_idx=inter;
            obj.outer_idx=indices;
            
        end
        
        function obj=set_merge_order(obj,merge_order)
            
            
            range=size(obj);
            dim=numel(range);
            if dim==2
                if find(range==1)
                    dim=1;
                end
                
            end
            if (~isnumeric(merge_order) || ~isvector(merge_order))&&~isempty(merge_order)
                error('Input numeric vector arrays or empty brackets to specify merge order');               
            end
            
            
            if sum(floor(merge_order)-merge_order)~=0 
                error('Merge indices must be integer valued')
            end
            if ~isempty(merge_order)&&(min(merge_order)<1 || max(merge_order)>dim)
                error(['Merger indices must be between 1 and ' int2str(dim)])
            end
           
            if isempty(merge_order)
                merge_order=1:dim;
            end
            
            indices=1:dim;            
            indices(merge_order)=0;            
            indices=find(indices);
            if indices
                error('Not all merger indices have been set');
            end
            obj.merge_idx=merge_order;
            
        end
        
        function obj=permute(obj,order)
            obj.data=permute(obj.data,order);
        
        end
        
        function sref = subsref(obj,s)
            switch s(1).type
                case '()'
                    sref=obj.data(s(1).subs{:});
                case '.'
                    sref = builtin('subsref',obj,s);
                otherwise
                    error('Invalid indexing');
                        
            end
        end
        
        function ret=sqrt(obj)
            ret=obj;
            ret.data=sqrt(ret.data);
        
        end
        
        function ret=uminus(obj)
            ret=obj;
            ret.data=-(ret.data);
        end
        
        A=flatten(A,row_idx,col_idx)
        C=mtimes(A,B)
        C=plus(A,B);
        C=minus(A,B);
        C=multinomial_solve(A,b)
        B=toLattice(A,row_idx,col_idx,depth_idx,dims);
    end
    
    methods (Access=private)
        C=back_solve(A)
    
    end
    
    methods (Access=protected)
        error_check_mult(A,B);
        error_check_mldivide(A,C);
        error_check_merge(A,B);
        error_check_flatten(A,row_idx,col_idx);
    end
end