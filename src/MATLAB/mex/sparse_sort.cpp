#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdint>
#include "SparseSortMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "LibMIAAlgorithm.h"
#include "check_sparse.h"


typedef long long index_type;

template<class T>
void perform_sort(T*a_data, index_type  * a_indices, mwSize a_size){
	using namespace LibMIA;
	internal::Introsort(a_indices, a_indices + a_size, std::less<index_type>(),
		internal::DualSwapper<index_type*, T*>(a_indices, a_data));

}

//!Assumes both are sorted
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	if (nrhs != 2)
		mexErrMsgTxt("Two arguments must be provided to sort function.");
	
	mwSize a_data_length;
	mxClassID a_id=check_sparse_unary(nrhs, prhs, &a_data_length);    
	
    
    switch (a_id)
    {
    case mxDOUBLE_CLASS:

        perform_sort((double *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]),a_data_length);
        break;
    case mxSINGLE_CLASS:
		perform_sort((float *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), a_data_length);
        break;
    case mxINT32_CLASS:
		perform_sort((int32_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), a_data_length);		
        break;
	case mxINT64_CLASS:
		perform_sort((int64_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), a_data_length);
		
		break;
	case mxUINT32_CLASS:
		perform_sort((uint32_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), a_data_length);
		
		break;
	case mxUINT64_CLASS:
		perform_sort((uint64_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), a_data_length);
		
		break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



}

