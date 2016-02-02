#include "SparseDenseLatticeMultMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"

#include "MappedDenseLattice.h"
#include "MappedSparseLattice.h"
#include "DenseLattice.h"
#include "LibMIAException.h"




template<class T>
int perform_solve(T * C, T*a_data, mex_index_type* a_index, const mex_index_type*a_subs, mwSize a_size, T*B, const mwSize*b_subs, const mwSize*c_subs){


	LibMIA::MappedSparseLattice<T> latA(a_data, a_index, a_size, a_subs[0], a_subs[1], a_subs[2]);
	const LibMIA::MappedDenseLattice<T> latB(B, b_subs[0], b_subs[1], b_subs[2]);	
	LibMIA::MappedDenseLattice<T> latC(C, c_subs[0], c_subs[1], c_subs[2]);
	
	int ret_code;

    try{
		
		latC = latA.solve(latB);
		//c_data = (T *)mxCalloc(latC.size(), sizeof(T));
		if (latC.solveInfo() == LibMIA::RankDeficient){
			ret_code = 0;
		}
		else if (latC.solveInfo() == LibMIA::FullyRanked){
			ret_code = 1;
		}
		else if (latC.solveInfo() == LibMIA::LeastSquares){
			ret_code = 2;
		}
		
    }
    catch(LibMIA::LatticeException& e){
        mexErrMsgTxt(e.what());
		return 0;
    }
	return ret_code;

}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	using namespace LibMIA::internal;
	if (nrhs != 5)
		mexErrMsgTxt("Four arguments must be provided to sparse solve dense.");
	if (nlhs != 1 && nlhs!=2)
		mexErrMsgTxt("One or two outputs required to solve a SparseLattice solve DenseLattice.");
	int dense_first = (int)check_passed_double(nrhs, prhs, 4);
	mwSize a_data_length,a_degree,a_index_length;
	mxClassID a_id,b_id;
	mwSize b_subs[3];
	mex_index_type *a_subs;
	mwSize c_subs[3];
	a_id = check_sparse_data(nrhs, prhs, 0, &a_data_length);
	a_index_length = check_sparse_indices(nrhs, prhs, 1);
	a_degree = check_sparse_dims(nrhs, prhs, 2);
	b_id = check_dense_params_lattice(nrhs, prhs, b_subs, 3);
	a_subs = (mex_index_type *)mxGetData(prhs[2]);
	if (!dense_first){
		
		if (b_subs[0] != a_subs[0]){
			mexErrMsgTxt("Left operand height must match right operand height.");
		}
		
		c_subs[0] = a_subs[1];
		c_subs[1] = b_subs[1];
		c_subs[2] = b_subs[2];
	}
	else{
		
		mexErrMsgTxt("Only function supported right now is sparse lattice inverted with dense lattice RHS.");
		/*if (a_subs[1] != b_subs[0]){
			mexErrMsgTxt("Left operand width must match right operand height.");
		}
		c_subs[0] = a_subs[0];
		c_subs[1] = b_subs[1];
		c_subs[2] = b_subs[2];*/

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
	if (a_subs[2] != b_subs[2]){
		mexErrMsgTxt("Depths of two lattices must be the same.");
	}
	

	 

	plhs[0] = mxCreateNumericArray(3, c_subs, a_id, mxREAL);
    
	int result;
    switch (a_id)
    {
    case mxDOUBLE_CLASS:
		
		result=perform_solve((double *)mxGetData(plhs[0]), (double *)mxGetData(prhs[0]), (mex_index_type*)mxGetData(prhs[1]), a_subs, a_data_length, (double *)mxGetData(prhs[3]), b_subs, c_subs);
        break;
    case mxSINGLE_CLASS:
		result = perform_solve((float *)mxGetData(plhs[0]), (float *)mxGetData(prhs[0]), (mex_index_type*)mxGetData(prhs[1]), a_subs, a_data_length, (float *)mxGetData(prhs[3]), b_subs, c_subs);
        break;
    
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized lattice data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }
	double*y;
	if (nlhs == 2){
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		y = mxGetPr(plhs[1]);
		*y = result;
	}


    //int fields =mxGetNumberOfFields(prhs[0]);
    //mexPrintf("Number of fields %d", fields);

}
