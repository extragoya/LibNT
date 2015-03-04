classdef (InferiorClasses = {?MIA}) SparseMIA <MIA
    
    properties
        
        mDims
        indices @int64 vector %linear indices, can be in arbitrary lexicographical order
        lexOrder
        isSorted
        
    end
    
    methods
        
        function obj=SparseMIA(varargin)
            
            if nargin ==0
                obj.mDims=[];
                obj.indices=int64([]);
                obj.data=[];
                obj.lexOrder=[];
                obj.isSorted=[];
                
            else
                if nargin==1
                    
                    arg=varargin{1};
                    if issparse(arg) %S = SparseLattice(data)
                        obj.mDims=size(arg);
                        obj.indices=int64(find(arg));
                        obj.data=nonzeros(arg);
                        obj.lexOrder=[1 2]; %equivalent to column major ordering
                        obj.isSorted=true;
                        
                    elseif isa(arg,'MIA') %S = SparseLattice(denseMIA)
                        obj.mDims=size(arg);
                        idx=find(arg.data);
                        obj.data=arg.data(idx);
                        obj.indices=int64(idx);
                        obj.lexOrder=1:numel(obj.dims);
                        obj.isSorted=true;                        
                    elseif ismatrix(arg) %assuming just a set of indices
                        if(size(arg,1)~=1)
                            error('Must input row vector for dimensions')
                        end
                        obj.mDims=arg;
                        obj.indices=int64([]);
                        obj.data=[];
                        obj.lexOrder=1:numel(obj.dims);
                        obj.isSorted=[];
                    else
                        error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    end
                elseif (nargin ==3 || nargin==4 || nargin==5)
                    %S = SparseLattice(data,indices,dims,lexOrder,isSorted)
                    obj.data=varargin{1};
                    obj.indices=int64(varargin{2});
                    if size(obj.indices,2)>1
                        error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    end
                    obj.mDims=varargin{3};
                    
                    idx= obj.dims>1;
                    obj.mDims=obj.mDims(idx);     
                    
                    
                    
                    obj.lexOrder=1:numel(obj.dims);
                    if nargin<4
                        obj.lexOrder=1:numel(obj.mDims);
                    else
                        obj.lexOrder=varargin{4};
                    end
                    if (nargin<5)
                        obj.isSorted=false;     
                    else  % if sorted is not specified, we assume it's unsorted
                        obj.isSorted=varargin{5};
                    end                            
                    
                    
                                  
                else
                    error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    
                end
                
            end
        end
        
        function ret=size(obj)
            
            ret=obj.dims;
            
        end
        function ret=nnz(obj)
            ret=length(obj.data);
        end
        function m_order=order(obj)
            m_order=length(obj.dims);
        
            
        end

        function ret=dims(obj)
            ret=obj.mDims;
        end
        
        function obj=insert(obj,data,index)
            if(~ismatrix(index))
                error('Must Input scalar or vector for index values');
            end
            if(size(index,2)==1)
                idx=index;
            else
                mult=obj.dims;
                mult=mult(obj.lexOrder); %change the dimensions to the current lexOrder
                mult=cumprod(mult); %compute cumulative product
                mult=[1 mult(1:end-1)];
                index=index(:,obj.lexOrder); %change indices to current lexOrder                
                idx=(index-1)*mult'+1;%compress them to linearized index
            end
            if(any(idx>obj.dimensionality()))
                error('Index must not be larger than dimensionality');
            end
            %adds arguments to object's data and indices arrays
            obj.data=[obj.data;data];
            obj.indices=[obj.indices;idx];
            obj.isSorted=false;
            
        end
        
        
        A=sort(A);
        A=permute(A,newLinIdx);
        A=changeLexOrder(A,newLinIdx);
        
        matrix = flatten(obj,row_inds,col_inds);
        
        C=do_plus(A,B,permute_idx);
        isequal=eq(a,b);
        is_eq = fuzzy_eq(A,B,tol);
        B=toLattice(A,row_idx,col_idx,depth_idx);
        
    end
    
    methods (Access=private)
        
        indices=pull_index(A,idx);      
        C=do_merge(A,B,permute_idx,op);
        
        
        
    end
    methods (Access=protected)
        A = assign(A,otherMIAA,assign_order);  
        is_equal=compare(A,B,func);
    end
    methods(Access=private, Static)
        new=reorder_from(old,index_order);
    end
end