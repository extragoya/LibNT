#include <iostream>
#include <sstream>
#include <string>
#include "SparseLatticeSolveMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"
#include "DenseLattice.h"
#include "MappedSparseLattice.h"
#include "LibMIAException.h"




template<class T>
mwSize perform_solve(mxArray *plhs[], mxClassID a_id, T*a_data, mex_index_type  * a_indices, T*b_data, mex_index_type * b_indices, mex_index_type *a_subs, mwSize a_size, mex_index_type *b_subs, mwSize b_size){

    
	T * c_data;
	

	LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, a_subs[0], a_subs[1], a_subs[2]);
	LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, b_subs[0], b_subs[1], b_subs[2]);

	try{
		LibMIA::DenseLattice<T> latC = latA.solve(latB);

		c_data = (T *)mxCalloc(latC.size(), sizeof(T));		
		std::copy(latC.data_begin(), latC.data_end(), &c_data[0]);		
		plhs[0] = mxCreateNumericMatrix(0, 0, a_id, mxREAL);
		mxSetData(plhs[0], c_data);
		mxSetM(plhs[0], latC.size());
		mxSetN(plhs[0], 1);		
		return latC.size();

	}
	catch (LibMIA::LatticeException& e){
		mexErrMsgTxt(e.what());
		return 0;

	}

    



}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	mex_index_type *a_subs;
	mex_index_type *b_subs;

	mwSize a_data_length;
	mwSize b_data_length;
	if (nlhs != 1)
		mexErrMsgTxt("One output arguments, c_data, representing dense output data.");
	else if (nlhs == 0)
		mexErrMsgTxt("No output");
	mxClassID a_id = check_sparse_params_lattice(nrhs, prhs, &a_subs, &b_subs, &a_data_length, &b_data_length);

    switch (a_id)
    {
    case mxDOUBLE_CLASS:

		perform_solve(plhs, a_id, (double *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), (double *)mxGetData(prhs[3]), (mex_index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length);
        break;
    case mxSINGLE_CLASS:
		perform_solve(plhs, a_id, (float *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), (float *)mxGetData(prhs[3]), (mex_index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length);
        break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized lattice data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



}

