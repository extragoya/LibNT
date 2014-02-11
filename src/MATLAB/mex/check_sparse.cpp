

#include "check_sparse.h"

mxClassID check_sparse_params(int nrhs, const mxArray *prhs[],double  *a_subs,double  *b_subs,mwSize* a_data_length,mwSize* b_data_length){

    if (nrhs!=6)
        mexErrMsgTxt("Six arguments must be provided to sparse lattice multiplication.");
    
    bool _valid= mxIsNumeric(prhs[0])& mxIsInt64 (prhs[1]) & mxIsDouble (prhs[2])& mxIsNumeric (prhs[3])& mxIsInt64 (prhs[4])&  mxIsDouble(prhs[5]);

    if (!_valid)
        mexErrMsgTxt("Invalid class. Indices must be int64 and dims must be double.\n");

    mxClassID a_id =mxGetClassID(prhs[0]);
    mxClassID b_id =mxGetClassID(prhs[3]);
    if (a_id!=b_id)
        mexErrMsgTxt("Lattice data must have the same underlying data type");



    mwSize a_data_dims=mxGetNumberOfDimensions(prhs[0]);
    mwSize a_index_dims=mxGetNumberOfDimensions(prhs[1]);
    mwSize a_dims_dims=mxGetNumberOfDimensions(prhs[2]);
    mwSize b_data_dims=mxGetNumberOfDimensions(prhs[3]);
    mwSize b_index_dims=mxGetNumberOfDimensions(prhs[4]);
    mwSize b_dims_dims=mxGetNumberOfDimensions(prhs[5]);
	_valid = a_data_dims == 2 && a_index_dims == 2 && b_data_dims == 2 && b_index_dims == 2 && a_dims_dims == 2 && b_dims_dims == 2;



    const mwSize* a_data_subs=mxGetDimensions(prhs[0]);
    const mwSize* a_index_subs=mxGetDimensions(prhs[1]);
    const mwSize* a_dims=mxGetDimensions(prhs[2]);
    const mwSize* b_data_subs=mxGetDimensions(prhs[3]);
    const mwSize* b_index_subs=mxGetDimensions(prhs[4]);
    const mwSize* b_dims=mxGetDimensions(prhs[5]);
	_valid = _valid && a_data_subs[1] == 1 && a_index_subs[1] == 1 && b_data_subs[1] == 1 && b_index_subs[1] == 1 && a_dims[0] == 1 && b_dims[0] == 1 && a_dims[1] == 3 && b_dims[1] == 3;
    if (!_valid)
        mexErrMsgTxt("Input must be two vectors representing a_data a_indices and a vector of dimensions [row cols tabs]. This is followed by the same info for b]");

    _valid= (a_data_subs[0]==a_index_subs[0]) & (b_data_subs[0]==b_index_subs[0]);
    if (!_valid)
        mexErrMsgTxt("Data and index vectors must be the same length.");
    
    std::copy(mxGetPr(prhs[2]),mxGetPr(prhs[2])+3,a_subs);
    std::copy(mxGetPr(prhs[5]),mxGetPr(prhs[5])+3,b_subs);

    *a_data_length=a_data_subs[0];
    *b_data_length=b_data_subs[0];
    return a_id;

}
