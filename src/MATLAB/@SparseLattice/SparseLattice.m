classdef SparseLattice < Lattice
    %Sparse Lattice class
    %Developed by Adam Harrison
    %U of A Electronic Imaging Lab, 2010
    
    %   S = SparseLattice(X) converts a sparse or full lattice(or matrix) to sparse form by
    %   squeezing out any zero elements.
    %
    %   S = SparseLattice(i,j,k,vals,m,n,p) uses the rows of [i,j,k,a] to generate an
    %   m-by-n-by-p sparse lattice.
    %   The three integer index vectors, i,j, and k and the real entries
    %   vector, a, all have the same length, nnz, which is the number of
    %   nonzeros in the resulting sparse lattice S .  Any elements of a
    %   which have duplicate values of i and j and k are added together.
    %
    %   S = SparseLattice(inds,vals,m,n,p) uses the rows of [inds,a] to generate an
    %   m-by-n-by-p sparse lattice.
    %   This corresponds to the internal storage of SparseLattice. inds
    %   contains linearized indices of values in a. Indices are stored in
    %   column, then row-major ordering.
    
    %   There are several simplifications of this seven argument call.%
    %
    %   S = SparseLattice(i,j,k,vals) uses m = max(i) and n = max(j) and p = max(k).
    %
    %   S = SparseLattice(m,n,p) abbreviates SparseLattice([],[],m,n,p).  This
    %   generates the ultimate sparse lattice, an m-by-n-by-p all zero lattice.
    
    
    
    
    
    
    
    
    properties (SetAccess = protected)
        
        
        
        inds
        indices
        
        
        
        nnz %number of non zeros
        
        iscompressed=false;
        
    end
    methods
        
        function obj = SparseLattice(varargin)
            obj@Lattice()
            if nargin==0
                obj.vals=[];
                obj.inds=int64([]);     
                obj.m=0;
                obj.n=0;
                obj.p=0;
                obj.nnz=0;
            elseif nargin==1
                arg=varargin{1};
                if isa(arg,'SparseLattice');
                    obj=arg;
                elseif isnumeric(arg) && numel(size(arg))>1 && numel(size(arg))<4
                    %create a spares lattice based from a full matrix or lattice
                    if numel(size(arg))==2
                        [t_m t_n]=size(arg);
                        t_p=1;
                        
                    else
                        [t_m t_n t_p]=size(arg);
                    end
                    obj.inds=int64(find(arg));
                    t_vals=arg(obj.inds);
                    obj.m=t_m;
                    obj.n=t_n;
                    obj.p=t_p;
                    obj.vals=t_vals;                              
                    obj.nnz=length(t_vals);
                else
                    error('Only 2 or 3-indexed numeric arrays may be converted to a SparseLattice object')
                end
            elseif nargin==7
                %S = SparseLattice(i,j,k,vals,m,n,p)
                
                obj.m=varargin{5};
                obj.n=varargin{6};
                obj.p=varargin{7};
                
                SparseLattice.error_check_range(obj.m,obj.n,obj.p);
                
                
                
                t_vals=varargin{4};
                
                t_vals=SparseLattice.error_check_nonzeros(t_vals);
                
                
                obj.nnz=length(t_vals);
                
                i=varargin{1};
                j=varargin{2};
                k=varargin{3};
                t_inds = int64(sub2ind(obj, i, j, k));
                t_inds=SparseLattice.error_check_idx(t_inds,obj.nnz,obj.m,obj.n,obj.p);     
                obj.inds=t_inds;
                
                [obj.inds,idx] = sort(obj.inds);
                t_vals=t_vals(idx);
                
                
                
                obj.vals=t_vals;         
                
               
            elseif nargin==5
                %S = SparseLattice(inds,vals,m,n,p)
                obj.m=varargin{3};
                obj.n=varargin{4};
                obj.p=varargin{5};
                
                SparseLattice.error_check_range(obj.m,obj.n,obj.p);
                
                t_vals=varargin{2};
                obj.nnz=length(t_vals);
                t_vals=SparseLattice.error_check_nonzeros(t_vals);
                
                
                t_inds=int64(varargin{1}); 
                t_inds=SparseLattice.error_check_idx(t_inds,obj.nnz,obj.m,obj.n,obj.p); 
                [obj.inds,idx]=sort(t_inds);
                t_vals=t_vals(idx);                
                obj.vals=t_vals;     
                
            elseif nargin==3
                %S = SparseLattice(m,n,p)
                obj.m=varargin{1};
                obj.n=varargin{2};
                obj.p=varargin{3};
                obj.vals=[];
                obj.inds =int64([]);                
                obj.nnz=0;
                
            else
                error('SparseLattice constructor called with incomptabile parameters. Please view help.')
            end
        end % SparseLattice
        
        function ret=size(obj)
            
            ret=[obj.m obj.n obj.p];
            
        end
        
        function disp(obj)
            
            
            disp('                 i                 j                 k                 val')
            
            disp([obj.row() obj.col() obj.depth() obj.vals])
        end
        
        
        
        function ret = subsref(obj,s)
            % SUBSREF Implementing the following syntax:
            % obj(i,j,k)
            % obj.col
            % obj.row
            % obj.depth
            switch s(1).type
                
                case '.'
                    switch s(1).subs
                        case 'col'
                            if strcmp(s(2).type,'()')
                                if(isempty(s(2).subs))
                                    ret = col(obj);
                                else
                                    ret = col(obj,s(2).subs{1});
                                end
                            else
                                error('Syntax for column is A.col(idx)')
                            end
                        case 'row'
                            if strcmp(s(2).type,'()')
                                if(isempty(s(2).subs))
                                    ret = row(obj);
                                else
                                    ret = row(obj,s(2).subs{1});
                                end
                            else
                                error('Syntax for row is A.row(idx)')
                            end
                        case 'tab'
                            if strcmp(s(2).type,'()')
                                if(isempty(s(2).subs))
                                    ret = tab(obj);
                                else
                                    ret = tab(obj,s(2).subs{1});
                                end
                            else
                                error('Syntax for depth is A.depth(idx)')
                            end
                        otherwise
                            ret = obj.(s.subs);
                            
                    end
                case '()'    
                    if length(s(1).subs)==1
                        ret=obj.vals(s(1).subs{1});

                    end
                    
                                    
                otherwise
                    error('Specify value for x as obj(x)')
            end
        end % subsref
        
        C=plus(A,B)
        C=minus(A,B)
        C=mtimes(A,B)
        
        A=transpose(A);
        A=ctranspose(A);
    end %public methods
    
    methods (Access=protected)
        C=union(A,B,fhandle)
        C=intersect(A,B,fhandle)
        C=binary_setup(A,B,fhandle,f_intersect)
        [A,B,IsSparse] = chkbinary(A,B)
        
        
        function [i j k]=ind2sub(obj,idx)
            if nargin==2
                [i j k]=ind2sub([obj.m,obj.n,obj.p],idx);
            else
                [i j k]=ind2sub([obj.m,obj.n,obj.p],obj.inds);
            end
        end
        
        function idx=sub2ind(obj,i,j,k)
            idx=sub2ind([obj.m,obj.n,obj.p],i,j,k);
        end
        
        function i=row(obj,idx)
            if nargin==1
                idx=obj.inds;
            end
            idx=int64(idx-1);
            i=mod(idx,obj.m)+1;            
        end
        
        function j=col(obj,idx)
            if nargin==1
                idx=obj.inds;
            end
            idx=int64(idx-1);
            j=idivide(idx,obj.m);
            j=mod(j,obj.n)+1;   
        end
        
        function k=tab(obj,idx)
            if nargin==1
                idx=obj.inds;
            end
           idx=int64(idx-1);
           k=idivide(idx,(obj.m*obj.n))+1;
           
        end
        
        
        
    end %protected methods
    
    methods (Access=protected, Static=true)
        
        [c c_inds]=merge_list(a,a_inds,b,b_inds,m,n,p);
    
   
        
        function inds=error_check_idx(inds,nnz,m,n,p)
            if ~isvector(inds)
                error('Input vector for row, column, and depth indices')
            end
            if length(inds)~=nnz
                error('Length of row, column, and depth indices must equal length of non-zero values')
            end
            if size(inds,1)==1
                inds=inds';
            end
           
                  
            if max(inds)>m*n*p
                error('Indices exceed valid range')
            end
            
           
        end
        
        function error_check_range(m,n,p)
            if ~isscalar(m) || ~(isa(m,'integer') || (imag(m)==0 && mod(m,1)==0) ) || m<1
                error('Input positive scalar integer value for number of rows')
            end
            
            if ~isscalar(n) || ~(isa(n,'integer') || (imag(n)==0 && mod(n,1)==0) ) || n<1
                error('Input positive scalar integer value for number of columns')
            end
            
            if ~isscalar(p) || ~(isa(p,'integer') || (imag(p)==0 && mod(p,1)==0) ) || p<1
                error('Input positive scalar integer value for depth size')
            end
            
        end
        
        function vals=error_check_nonzeros(vals)
            if ~isvector(vals)
                error('Input vector for nonzero values')
            end
            if size(vals,1)==1
                vals=vals';
            end
        end
        
    end %protected, static methods 
    
end % classdef
