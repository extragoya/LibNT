#include "SparseDenseLatticeMultMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"

#include "MappedDenseLattice.h"
#include "MappedSparseLattice.h"
#include "DenseLattice.h"
#include "LibMIAException.h"




template<class T>
mwSize perform_mult(mxArray *plhs[], mxClassID a_id, T*a_data, mex_index_type* a_index, const mex_index_type*a_subs, mwSize a_size, T*B, const mwSize*b_subs, int dense_first){


	LibMIA::MappedSparseLattice<T> latA(a_data, a_index, a_size, a_subs[0], a_subs[1], a_subs[2]);
	const LibMIA::MappedDenseLattice<T> latB(B, b_subs[0], b_subs[1], b_subs[2]);
	LibMIA::SparseLattice<T> latC;
	T * c_data;
	mex_index_type * c_indices;

    try{
		if (dense_first)
			latC = latB*latA;
		else
			latC = latA*latB;
		c_data = (T *)mxCalloc(latC.size(), sizeof(T));
		c_indices = (mex_index_type *)mxCalloc(latC.size(), sizeof(mex_index_type));

		std::copy(latC.data_begin(), latC.data_end(), &c_data[0]);
		std::copy(latC.index_begin(), latC.index_end(), &c_indices[0]);
		plhs[0] = mxCreateNumericMatrix(0, 0, a_id, mxREAL);
		mxSetData(plhs[0], c_data);
		mxSetM(plhs[0], latC.size());
		mxSetN(plhs[0], 1);
		plhs[1] = mxCreateNumericMatrix(0, 0, mxINT64_CLASS, mxREAL);
		mxSetData(plhs[1], c_indices);
		mxSetM(plhs[1], latC.size());
		mxSetN(plhs[1], 1);
		return latC.size();
		
    }
    catch(LibMIA::LatticeException& e){
        mexErrMsgTxt(e.what());
		return 0;
    }

}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	using namespace LibMIA::internal;
	if (nrhs != 5)
		mexErrMsgTxt("Five arguments must be provided to sparse/dense.");
	if (nlhs != 2)
		mexErrMsgTxt("Two output arguments, c_data c_indices, required to multiply Dense and SparseLattices.");
	int dense_first = (int)check_passed_double(nrhs, prhs, 4);
	mwSize a_data_length,a_degree,a_index_length;
	mxClassID a_id,b_id;
	mwSize b_subs[3];
	mex_index_type *a_subs;
	
	a_id = check_sparse_data(nrhs, prhs, 0, &a_data_length);
	a_index_length = check_sparse_indices(nrhs, prhs, 1);
	a_degree = check_sparse_dims(nrhs, prhs, 2);
	b_id = check_dense_params_lattice(nrhs, prhs, b_subs, 3);
	a_subs = (mex_index_type *)mxGetData(prhs[2]);
	if (dense_first){
		
		if (b_subs[1] != (mwSize)a_subs[0]){
			mexErrMsgTxt("Left operand width must match right operand height.");
		}
		
	}
	else{
		
		
		if ((mwSize)a_subs[1] != b_subs[0]){
			mexErrMsgTxt("Left operand width must match right operand height.");
		}
		

	}
	if (a_degree != 3){
		mexErrMsgTxt("Three dimensions must be passed in for sparse dimensions.");
	}
	if (a_index_length != a_data_length){
		mexErrMsgTxt("Sparse index vector and data vector must be the same length.");
	}
	if (a_id != b_id){
		mexErrMsgTxt("Operands must be of the same datatype.");
	}
	if ( (mwSize)a_subs[2] !=b_subs[2]){
		mexErrMsgTxt("Depths of two lattices must be the same.");
	}
	

	 

   
    //plhs[0]=mxCreateNumericArray(0, 0,   a_id, mxREAL);

    switch (a_id)
    {
    case mxDOUBLE_CLASS:
		
		perform_mult(plhs, a_id,(double *)mxGetData(prhs[0]), (mex_index_type*)mxGetData(prhs[1]), a_subs, a_data_length, (double *)mxGetData(prhs[3]), b_subs, dense_first);
        break;
    case mxSINGLE_CLASS:
		perform_mult(plhs, a_id, (float *)mxGetData(prhs[0]), (mex_index_type*)mxGetData(prhs[1]), a_subs, a_data_length, (float *)mxGetData(prhs[3]), b_subs, dense_first);
        break;
    case mxINT32_CLASS:
		perform_mult(plhs, a_id, (int32_t *)mxGetData(prhs[0]), (mex_index_type*)mxGetData(prhs[1]), a_subs, a_data_length, (int32_t *)mxGetData(prhs[3]), b_subs, dense_first);
        break;
	case mxINT64_CLASS:
		perform_mult(plhs, a_id, (int64_t *)mxGetData(prhs[0]), (mex_index_type*)mxGetData(prhs[1]), a_subs, a_data_length, (int64_t *)mxGetData(prhs[3]), b_subs, dense_first);
		break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized lattice data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



    //int fields =mxGetNumberOfFields(prhs[0]);
    //mexPrintf("Number of fields %d", fields);

}
