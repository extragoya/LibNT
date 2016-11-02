

#include "check_sparse.h"



mxClassID check_sparse_params_lattice(int nrhs, const mxArray *prhs[], mex_index_type  *a_subs[], mex_index_type  *b_subs[], mwSize* a_data_length, mwSize* b_data_length,bool allowed_seven){

	if (!allowed_seven){
		if (nrhs != 6)
			mexErrMsgTxt("Wrong number of arguments. Input must be two vectors representing a_data a_indices and a vector of dimensions [row cols tabs]. This is followed by the same info for b");
	}
	else{
		if (nrhs < 6 && nrhs>7)
			mexErrMsgTxt("Wrong number of arguments. Input must be two vectors representing a_data a_indices and a vector of dimensions [row cols tabs]. This is followed by the same info for b. You are optionally allowed a seventh argrument indicating the multiplication to perform");
	}
    
    
	mwSize a_index_length,b_index_length;
	mxClassID a_id=check_sparse_data(nrhs, prhs, 0, a_data_length);
	a_index_length=check_sparse_indices(nrhs, prhs, 1);
	mxClassID b_id = check_sparse_data(nrhs, prhs, 3, b_data_length);
	b_index_length = check_sparse_indices(nrhs, prhs, 4);
	mwSize a_degree = check_sparse_dims(nrhs, prhs, 2);
	mwSize b_degree = check_sparse_dims(nrhs, prhs, 5);
	
	bool _valid = a_degree == 3 && b_degree == 3;
    if (!_valid)
        mexErrMsgTxt("Must input height, width, and depth dimensions for each input.\n");
    
    
    if (a_id!=b_id)
        mexErrMsgTxt("Lattice data must have the same underlying data type");
  
        

	_valid = (a_index_length == *a_data_length) && (b_index_length == *b_data_length);
    if (!_valid)
        mexErrMsgTxt("Data and index vectors must be the same length.");
    
	*a_subs = (mex_index_type *)mxGetData(prhs[2]);
	*b_subs = (mex_index_type *)mxGetData(prhs[5]);
    
    return a_id;

}

mxClassID check_dense_params_lattice(int nrhs, const mxArray *prhs[], mwSize  *a_subs, int param_index){

	assert(nrhs > param_index);

	bool _valid = mxIsNumeric(prhs[param_index]);

	if (!_valid)
		mexErrMsgTxt("Invalid class\n");

	mxClassID a_id = mxGetClassID(prhs[param_index]);
	
	mwSize a_subs_l = mxGetNumberOfDimensions(prhs[param_index]);
	
	if (a_subs_l != 3 && a_subs_l != 2)
		mexErrMsgTxt("Lattice data must be three-dimensional");

	const mwSize* a_subs_m = mxGetDimensions(prhs[param_index]);
	

	
	a_subs[0] = a_subs_m[0];
	a_subs[1] = a_subs_m[1];
	a_subs[2] = 1;

	
	if (a_subs_l == 3){
		a_subs[2] = a_subs_m[2];		
	}

	return a_id;

}

double check_passed_double(int nrhs, const mxArray *prhs[], int param_index)
{
	assert(nrhs > param_index);

	bool _valid = mxIsDouble(prhs[param_index]);

	if (!_valid)
		mexErrMsgIdAndTxt("MEX:incorrectArguments","Expecting double in parameter %d.\n", param_index);

	mwSize op_dims = mxGetNumberOfDimensions(prhs[param_index]);

	_valid=op_dims == 2;
	if (!_valid)
		mexErrMsgIdAndTxt("MEX:incorrectArguments", "Expecting scalar double in parameter %d.\n", param_index);
	_valid = mxGetM(prhs[param_index]) == 1 && mxGetN(prhs[param_index]) == 1;
	if (!_valid)
		mexErrMsgIdAndTxt("MEX:incorrectArguments", "Expecting scalar double in parameter %d.\n", param_index);

	return mxGetScalar(prhs[param_index]);

}



mxClassID check_sparse_params_merge(int nrhs, const mxArray *prhs[],  mwSize* a_data_length, mwSize* b_data_length,double * op){

	if (nrhs != 5)
		mexErrMsgTxt("Five arguments must be provided to sparse lattice merge.");

	mwSize a_index_length, b_index_length;
	mxClassID a_id = check_sparse_data(nrhs, prhs, 0, a_data_length);
	a_index_length = check_sparse_indices(nrhs, prhs, 1);
	check_sparse_data(nrhs, prhs, 2, b_data_length);
	b_index_length = check_sparse_indices(nrhs, prhs, 3);
	check_passed_double(nrhs, prhs, 4);
	bool _valid = (a_index_length == *a_data_length) && (b_index_length == *b_data_length);
	if (!_valid)
		mexErrMsgTxt("Data and index vectors must be the same length.");
	
	
	*op = std::floor(mxGetScalar(prhs[4]));
	if (*op<0 || *op>1)
		mexErrMsgTxt("Input one or two for operation, ie 0 for add, 1 for subtract");

	
	return a_id;

}

mxClassID check_sparse_unary(int nrhs, const mxArray *prhs[], mwSize* a_data_length)
{
	assert(nrhs > 1);

	bool _valid = mxIsNumeric(prhs[0])& mxIsInt64(prhs[1]);

	if (!_valid)
		mexErrMsgTxt("Invalid class. Indices must be int64 and dims must be double.\n");

	mxClassID a_id = mxGetClassID(prhs[0]);



	mwSize a_data_dims = mxGetNumberOfDimensions(prhs[0]);
	mwSize a_index_dims = mxGetNumberOfDimensions(prhs[1]);
	
	_valid = a_data_dims == 2 && a_index_dims == 2;



	const mwSize* a_data_subs = mxGetDimensions(prhs[0]);
	const mwSize* a_index_subs = mxGetDimensions(prhs[1]);	

	_valid = _valid && (a_data_subs[1] == 1 || a_data_subs[0] == 0) && (a_index_subs[1] == 1 || a_index_subs[0] == 0);
	if (!_valid)
		mexErrMsgTxt("Input must be two vectors representing data and then indices.");

	_valid = (a_data_subs[0] == a_index_subs[0]);
	if (!_valid)
		mexErrMsgTxt("Data and index vectors must be the same length.");


	*a_data_length = a_data_subs[0];
	
	return a_id;
}

mwSize check_sparse_indices(int nrhs, const mxArray *prhs[],int param_index)
{
	assert(nrhs > param_index);

	bool _valid = mxIsInt64(prhs[param_index]);

	if (!_valid)
		mexErrMsgTxt("Invalid class. Indices must be int64.\n");	
	
	mwSize a_index_dims = mxGetNumberOfDimensions(prhs[param_index]);
	_valid = a_index_dims == 2;
	
	const mwSize* a_index_subs = mxGetDimensions(prhs[param_index]);

	_valid = _valid && (a_index_subs[1] == 1 || a_index_subs[0] == 0);
	if (!_valid)
		mexErrMsgTxt("Input must a vector representing indices.");

	return a_index_subs[0];

	
}

mxClassID check_sparse_data(int nrhs, const mxArray *prhs[], int param_index, mwSize* data_length)
{
	assert(nrhs > param_index);

	mxClassID a_id = mxGetClassID(prhs[param_index]);
	bool _valid = mxIsNumeric(prhs[param_index]);
	if (!_valid)
		mexErrMsgTxt("Invalid class. Data must be numeric.\n");

	mwSize a_data_dims = mxGetNumberOfDimensions(prhs[param_index]);
	_valid = a_data_dims == 2;

	const mwSize* a_data_subs = mxGetDimensions(prhs[param_index]);

	_valid = _valid && (a_data_subs[1] == 1 || a_data_subs[0] == 0);
	if (!_valid)
		mexErrMsgTxt("Input must a vector representing data.");
	*data_length = a_data_subs[0];
	return a_id;


}

mwSize check_sparse_linIdx(int nrhs, const mxArray *prhs[], int param_index )
{
	assert(nrhs > param_index+1); //must be enought inputs for current_linIdx desired_linIdx and dimensions

	bool _valid = mxIsNumeric(prhs[param_index]) && mxIsNumeric(prhs[param_index+1]);

	if (!_valid)
		mexErrMsgTxt("Invalid class. Permute indices must be numeric and dimensions must be int64.\n");

	mwSize a_currentlinIdx_dims = mxGetNumberOfDimensions(prhs[param_index]);
	mwSize a_desiredlinIdx_dims = mxGetNumberOfDimensions(prhs[param_index+1]);
	
	_valid = a_currentlinIdx_dims == 2 && a_desiredlinIdx_dims == 2;

	const mwSize* a_currentlinIdx_subs = mxGetDimensions(prhs[param_index]);
	const mwSize* a_desiredlinIdx_subs = mxGetDimensions(prhs[param_index + 1]);

	_valid = _valid && (a_currentlinIdx_subs[0] == 1 || a_currentlinIdx_subs[1] == 0) && (a_desiredlinIdx_subs[0] == 1 || a_desiredlinIdx_subs[1] == 0);
	if (!_valid)
		mexErrMsgTxt("Input must row vectors representing current and desired lexicographical order.");
	_valid = _valid && a_currentlinIdx_subs[1] == a_desiredlinIdx_subs[1];
	if (!_valid)
		mexErrMsgTxt("Vectors of lexicographical order must be the same length.");

	mxClassID a_id = mxGetClassID(prhs[param_index]);
	mxClassID b_id = mxGetClassID(prhs[param_index+1]);
	_valid = _valid && a_id == mxUINT8_CLASS && b_id == mxUINT8_CLASS;
	if (!_valid)
		mexErrMsgTxt("Vectors of lexicographical order must be unsigned integer.");

	return a_currentlinIdx_subs[1];


}

mwSize check_sparse_dims(int nrhs, const mxArray *prhs[], int param_index){

	assert(nrhs > param_index);

	bool _valid = mxIsInt64(prhs[param_index]);

	if (!_valid)
		mexErrMsgTxt("Dimensions must be int64.\n");



	mwSize a_dim_dims = mxGetNumberOfDimensions(prhs[param_index]);
	

	_valid = a_dim_dims == 2;



	const mwSize* a_dim_subs = mxGetDimensions(prhs[param_index]);
	

	_valid = _valid && (a_dim_subs[0] == 1 || a_dim_subs[1] == 0);
	if (!_valid)
		mexErrMsgTxt("Input must be row vector of dimensions.");		

	return a_dim_subs[1];

}
