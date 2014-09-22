classdef (InferiorClasses = {?MIA}) SparseMIA <MIA
    
    properties
        
        dims
        indices %linear indices, can be in arbitrary lexicographical order
        linIdx
        isSorted
        
    end
    
    methods
        
        function obj=SparseMIA(varargin)
            
            if nargin ==0
                obj.dims=[];
                obj.indices=[];
                obj.data=[];
                obj.linIdx=[];
                obj.isSorted=[];
                
            else
                if nargin==1
                    arg=varargin{1};
                    if issparse(arg)
                        obj.dims=size(arg);
                        obj.indices=int64(find(arg));
                        obj.data=nonzeros(arg);
                        obj.linIdx=[1 2]; %equivalent to column major ordering
                        obj.isSorted=true;
                        
                    elseif isa(arg,'MIA')
                        obj.dims=size(arg);
                        idx=find(arg.data);
                        obj.data=arg.data(idx);
                        obj.indices=int64(idx);
                        obj.linIdx=1:numel(obj.dims);
                        obj.isSorted=true;                        
                        
                    else
                        error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    end
                elseif (nargin ==3 || nargin==4)
                    %TODO error checking
                    obj.data=varargin{1};
                    obj.indices=int64(varargin{2});
                    if size(obj.indices,2)>1
                        error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    end
                    obj.dims=varargin{3};
                    
                    idx= obj.dims>1;
                    obj.dims=obj.dims(idx);     
                    
                    
                    
                    obj.linIdx=1:numel(obj.dims);
                    if nargin==3 % if sorted is not specified, we assume it's unsorted
                        obj.isSorted=false;     
                    else
                        obj.isSorted=varargin{4};
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
        %TODO mldivide
        A=sort(A,newLinIdx);
        B=flatten(A,row_idx,col_idx);
        C=mtimes(A,B)
        C=plus(A,B)
        [indices parition]=extract_indices(A,idx,linear);
        indices=pull_index(A,idx);
        A=consolidate_indices(A);
    end
    
    methods (Access=private)
        C=perform_merge(A,B,fhandle);
        C=union(A,B,fhandle);
        
    end
end