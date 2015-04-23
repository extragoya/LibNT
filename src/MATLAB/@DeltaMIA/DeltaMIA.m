classdef (InferiorClasses = {?MIA,?SparseMIA}) DeltaMIA <SparseMIA
    
    methods
        
        function obj=DeltaMIA(varargin)
            doerror=false;
            if nargin==1
                arg=varargin{1};
                if isscalar(arg)
                    if floor(arg)==arg
                        obj.dims=[arg arg];
                        indices=1:arg;
                        obj.indices=((indices-1)*arg+indices)';
                        obj.inter_idx=[1 2];
                        obj.inner_idx=[];
                        obj.outer_idx=[];
                        obj.data=ones(arg,1);
                        obj.partition=numel(obj.dims);
                    else
                        doerror=true;
                    end
                else
                    doerror=true;
                end
                
            else
                doerror=true;
            end
            if doerror
                error('DeltaMIA constructor called with invalid argument, see help');
            end
        end
        C=mtimes(A,B)
    end
end