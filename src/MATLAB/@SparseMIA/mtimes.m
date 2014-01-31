function C=mtimes(A,B)
%SparseMIA multiplication
%
%C=A*B multiplies two SparseMIAs together. Indices undergoing inner, inter,
%and outer product must be set by member function set_idx. If indices are
%not set, A and B undergo by default an inter product for each
%corresponding element.
%
%Order of indices after multiplication are [A.outer_idx inter_idx B.outer_idx]
%
%Function operates by mapping A and B into lattices based on which indices
%are set for inner, inter, and outer product. To avoid lattice indices that
%potentially suffer from integer overflow, A and B are mapped into
%compressed lattices.
%
%Developed by Adam Harrison
%U of A Electronic Imaging Lab, 2011

if ~isa(A,'MIA') || ~isa(B,'MIA')
    error('SparseMIAs can only be multipled with classes or subclasses of MIA');
end

error_check_mult(A,B); %ensure indices match up properly

switched=false;
if ~isa(A,'SparseMIA')
    temp=A;
    A=B;
    B=temp;
    clear temp
    switched=true;
end

%extract multiplication indices

a_outer_idx=A.outer_idx;
b_outer_idx=B.outer_idx;
a_inter_idx=A.inter_idx;
b_inter_idx=B.inter_idx;
a_inner_idx=A.inner_idx;
b_inner_idx=B.inner_idx;

a_length=length(A.data);

%extract indice ranges for resulting MIA C
a_dims=A.size;
if numel(a_dims==2)
    if a_dims(2)==1
        a_dims=a_dims(1);
    end
end
b_dims=B.size;
if numel(b_dims==2)
    if b_dims(2)==1
        b_dims=b_dims(1);
    end
end
a_row_dims=a_dims(a_outer_idx);
a_col_dims=a_dims(a_inner_idx);
b_col_dims=b_dims(b_outer_idx);
b_tab_dims=b_dims(b_inter_idx);



%extact outer product indices of each element for A and calculate a
%one-to-one linear and compressed mapping
if ~isempty(a_outer_idx)
    [a_row_map a_row_partition]=A.extract_indices(a_outer_idx);
    [a_row_map,~,a_row_back_map]=unique(a_row_map,'rows');
    a_m=size(a_row_map,1); %number of rows in lattice version of A
else
    a_row_map=[];
    a_row_partition=[];
    a_row_back_map=ones(a_length,1);
    a_m=1;
end
a_row_back_map=int64(a_row_back_map);


if isa(B,'SparseMIA')
    %extact outer product indices of each element for B and calculate a
    %one-to-one linear and compressed linear mapping
    b_length=length(B.data);
    
    if ~isempty(b_outer_idx)
        [b_col_map b_col_partition]=B.extract_indices(b_outer_idx);
        [b_col_map,~,b_col_back_map]=unique(b_col_map,'rows');
        b_n=size(b_col_map,1); %number of columns in lattice version of B
    else
        
        b_col_map=[];
        b_col_partition=[];
        b_col_back_map=ones(b_length,1);
    end
    
    %extact inner product indices for each element of A and B. Calculate a
    %one-to-one linear and compressed linear mapping using both sets of indices
    if ~isempty(a_inner_idx)
        [a_inner_map inner_parition]=A.extract_indices(a_inner_idx);
        b_inner_map=B.extract_indices(b_inner_idx);
        inner_map=[a_inner_map;b_inner_map];
        
        [inner_map,~,inner_back_map]=unique(inner_map,'rows');
        a_col_back_map=inner_back_map(1:a_length);
        b_row_back_map=inner_back_map(a_length+1:end);
        
        a_n=size(inner_map,1);
    else
        a_col_back_map=ones(a_length,1);
        b_row_back_map=ones(b_length,1);
        inner_parition=[];
        a_n=1;
    end
    b_m=a_n; %number of columns and rows in lattice versions of A and B
    
    clear inner_map a_inner_map b_inner_map
    
    %perform identical operation on inter products as done for inner products
    if ~isempty(a_inter_idx)
        [a_inter_map inter_partition]=A.extract_indices(a_inter_idx);
        b_inter_map=B.extract_indices(b_inter_idx);
        inter_map=[a_inter_map;b_inter_map];
        
        [inter_map,~,inter_back_map]=unique(inter_map,'rows');
        a_tab_back_map=inter_back_map(1:a_length);
        b_tab_back_map=inter_back_map(a_length+1:end);
        
        a_p=size(inter_map,1);
    else
        a_tab_back_map=ones(a_length,1);
        b_tab_back_map=ones(b_length,1);
        inter_partition=[];
        a_p=1;
        inter_map=ones(1,1);
    end
    a_col_back_map=int64(a_col_back_map);
    a_row_back_map=int64(a_row_back_map);
    a_tab_back_map=int64(a_tab_back_map);
    b_row_back_map=int64(b_row_back_map);
    b_col_back_map=int64(b_col_back_map);
    b_tab_back_map=int64(b_tab_back_map);
    
    
    b_p=a_p;
    
    
    %flatten A and B into lattices and multiply
    A_lat=SparseLattice(a_row_back_map,a_col_back_map,a_tab_back_map,A.data,a_m,a_n,a_p);
    B_lat=SparseLattice(b_row_back_map,b_col_back_map,b_tab_back_map,B.data,b_m,b_n,b_p);
    
    
else %B is an MIA - no need to perform B outer product mappings or to find common inner product indices between A and B
    b_col_dims=b_dims(b_outer_idx);
    b_col_partition=numel(b_col_dims);
    if b_col_partition==0
        b_col_partition=[];
    end
    b_n=prod(b_col_dims);
    b_row_dims=b_dims(b_inner_idx);
    b_m=prod(b_row_dims);
    a_n=b_m;
    b_col_map=(1:b_n)';
    
    
    
    if ~isempty(a_inner_idx) %indices shouldn't exceed maxint since they share the same range as MIA B
        a_col_back_map=int64(A.extract_indices(a_inner_idx));   
        
    else
        a_col_back_map=int64(ones(a_length,1));
        
    end
    
    if ~isempty(a_inter_idx)
        [a_tab_back_map]=int64(A.extract_indices(a_inter_idx));
        inter_partition=numel(a_inter_idx);
        b_tab_dims=b_dims(b_inter_idx);
        b_p=prod(b_tab_dims);
        
        inter_map=(1:b_p)';
    else
        a_tab_back_map=int64(ones(a_length,1));
        b_p=1;
        inter_map=[];
        inter_partition=[];
    end
    a_p=b_p;
    
    
    A_lat=SparseLattice(a_row_back_map,a_col_back_map,a_tab_back_map,A.data,a_m,a_n,a_p);
    B_lat=B.toLattice(b_inner_idx, b_outer_idx, b_inter_idx,[b_m b_n b_p]);
    
    
end
clear a_inter_map b_inter_map

C_lat=A_lat*B_lat;
clear A_lat B_lat



c_data=C_lat.vals;
c_indices=[];
c_partition=[];

%calculate uncompressed MIA indices of C and convert to Sparse MIA

if switched
    if(~isempty(b_outer_idx))
        c_cols=C_lat.col();
        c_cols=b_col_map(c_cols);
        c_indices=[c_indices c_cols];
        c_partition=[c_partition length(b_outer_idx)];
    end
    
    if(~isempty(a_inter_idx))
        c_tabs=C_lat.tab();
        c_tabs=inter_map(c_tabs);
        c_indices=[c_indices c_tabs];
        c_partition=[c_partition length(a_inter_idx)];
    end
    
    if(~isempty(a_outer_idx))
        c_rows=C_lat.row();
        c_rows=a_row_map(c_rows);
        c_indices=[c_indices c_rows];
        c_partition=[c_partition length(a_outer_idx)];
    end
    for i=2:length(c_partition)
        c_partition(i)=c_partition(i)+c_partition(i-1);
    end    
    

    c_dims=[b_col_dims   b_tab_dims  a_row_dims];
    
   
    C=SparseMIA(c_data,c_indices,c_dims, c_partition);
    
else
    
    if(~isempty(a_outer_idx))
        c_rows=C_lat.row();
        c_rows=a_row_map(c_rows);
        c_indices=[c_indices c_rows];
        c_partition=[c_partition length(a_outer_idx)];
    end
    
    
    
    if(~isempty(a_inter_idx))
        c_tabs=C_lat.tab();
        c_tabs=inter_map(c_tabs);
        c_indices=[c_indices c_tabs];
        c_partition=[c_partition length(a_inter_idx)];
    end
    
    if(~isempty(b_outer_idx))
        c_cols=C_lat.col();
        c_cols=b_col_map(c_cols);
        c_indices=[c_indices c_cols];
        c_partition=[c_partition length(b_outer_idx)];
    end
    
    for i=2:length(c_partition)
        c_partition(i)=c_partition(i)+c_partition(i-1);
    end    
    

    c_dims=[a_row_dims b_tab_dims b_col_dims     ];
    
   
    C=SparseMIA(c_data,c_indices,c_dims, c_partition);
end
C=C.consolidate_indices;

