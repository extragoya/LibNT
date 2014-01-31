classdef (InferiorClasses = {?MIA}) SparseMIA <MIA
    
    properties
        
        dims
        indices %can be two dimensional based on partition array
        partition
        
    end
    
    methods
        
        function obj=SparseMIA(varargin)
            
            if nargin ==0
                obj.dims=[];
                obj.indices=[];
                obj.data=[];
                obj.inner_idx=[];
                obj.outer_idx=[];
                obj.inter_idx=[];
            else
                if nargin==1
                    arg=varargin{1};
                    if issparse(arg)
                        obj.dims=size(arg);
                        obj.indices=find(arg);
                        obj.data=nonzeros(arg);
                        obj.partition=2;
                        obj.inter_idx=1:2;
                        obj.merge_idx=1:2;
                        obj.inner_idx=[];
                        obj.outer_idx=[];
                    elseif isa(arg,'MIA')
                        obj.dims=size(arg);
                        idx=find(arg.data);
                        obj.data=arg.data(idx);
                        obj.indices=idx;
                        obj.partition=numel(obj.dims);
                        
                        if obj.partition==2
                            if obj.dims(2)==1
                                obj.partition=1;
                                
                            end
                        end
                        obj.inter_idx=arg.inter_idx;
                        obj.inner_idx=arg.inner_idx;
                        obj.outer_idx=arg.outer_idx;
                        obj.merge_idx=arg.merge_idx;
                    else
                        error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    end
                elseif nargin ==3
                    %TODO error checking
                    obj.data=varargin{1};
                    obj.indices=int64(varargin{2});
                    if size(obj.indices,2)>1
                        error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    end
                    obj.dims=varargin{3};
                    
                    idx= obj.dims>1;
                    obj.dims=obj.dims(idx);
                    l=numel(obj.dims);
                    if numel(obj.dims)==1
                        obj.dims=[obj.dims 1];
                        l=1;
                    end
                    obj.partition=numel(obj.dims);
                    [obj.indices idx] = sort(obj.indices);
                    obj.data=obj.data(idx);
                    obj.inter_idx=1:l;
                    obj.merge_idx=1:l;
                    obj.inner_idx=[];
                    obj.outer_idx=[];
                    
                elseif nargin ==4
                    %TODO error checking and cleaning up dimensions of unit
                    %size
                    obj.data=varargin{1};
                    obj.indices=int64(varargin{2});
                    obj.dims=varargin{3};
                    obj.partition=varargin{4};
                    
                    
                    [obj.indices idx] = sortrows(obj.indices);
                    obj.data=obj.data(idx);
                    l=numel(obj.dims);
                    obj.inter_idx=1:l;
                    obj.merge_idx=1:l;
                    obj.inner_idx=[];
                    obj.outer_idx=[];
                else
                    error('SparseMIA constructor called with incomptabile parameters. Please view help.')
                    
                end
                
            end
        end
        
        function ret=size(obj)
            
            ret=obj.dims;
            
        end
        %TODO mldivide
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