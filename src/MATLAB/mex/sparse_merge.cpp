#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "SparseMergeMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"
#include "SparseLattice.h"
#include "MappedSparseLattice.h"
#include "LibMIAException.h"


typedef long long index_type;

template<class T>
mwSize perform_merge(mxArray *plhs[],mxClassID a_id, T*a_data, index_type  * a_indices,  T*b_data, index_type * b_indices, mwSize a_size,mwSize b_size,double op){

    
	
	T * c_data;
    index_type  * c_indices;
	//create empty matrices
	plhs[0] = mxCreateNumericMatrix(0, 0, a_id, mxREAL);
	plhs[1] = mxCreateNumericMatrix(0, 0, mxINT64_CLASS, mxREAL);
	mwSize c_size=0;
	if (!a_size && !b_size){
		return 0;
	}
	else if (!b_size){
		c_data = (T *)mxCalloc(a_size, sizeof(T));
		c_indices = (index_type *)mxCalloc(a_size, sizeof(index_type));
		if (!c_data || !c_indices){
			mexErrMsgTxt("Not enough memory to allocate merge result");
			return 0;
		}
		std::copy(&a_data[0], &a_data[0]+a_size, &c_data[0]);
		std::copy(&a_indices[0], &a_indices[0] + a_size, &c_indices[0]);
		c_size= a_size;
	}
	else if (!a_size){
		c_data = (T *)mxCalloc(b_size, sizeof(T));
		c_indices = (index_type *)mxCalloc(b_size, sizeof(index_type));
		if (!c_data || !c_indices){
			mexErrMsgTxt("Not enough memory to allocate merge result");
			return 0;
		}
		if (op == 0){
			std::copy(&b_data[0], &b_data[0] + b_size, &c_data[0]);
		}
		else{
			auto it2 = c_data;
			for (auto it = b_data; it < b_data + b_size;++it2,++it){
				*it2 = -1*(*it);
			}
		}
		std::copy(&b_indices[0], &b_indices[0] + b_size, &c_indices[0]);
		c_size= b_size;
	}
	else{

		
		mwSize dummy_size = std::max(*(a_indices + a_size - 1), *(b_indices + b_size - 1));
		LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, dummy_size, 1, 1); //dummy lattices
		LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, dummy_size, 1, 1);
		LibMIA::SparseLattice<T> latC;






		try{

			if (op == 0)
				latC = latA + latB;
			else{
				latC = latA - latB;
			}
			c_size = latC.size();
			c_data = (T *)mxCalloc(c_size, sizeof(T));
			c_indices = (index_type *)mxCalloc(c_size, sizeof(index_type));
			if (!c_data || !c_indices){
				mexErrMsgTxt("Not enough memory to allocate merge result");
				return 0;
			}
			std::copy(latC.data_begin(), latC.data_end(), &c_data[0]);
			std::copy(latC.index_begin(), latC.index_end(), &c_indices[0]);

			
		}
		catch (LibMIA::LatticeException& e){
			mexErrMsgTxt(e.what());
			return 0;

		}
	}
    
	mxSetData(plhs[0], c_data);
	mxSetM(plhs[0], c_size);
	mxSetN(plhs[0], 1);
	mxSetData(plhs[1], c_indices);
	mxSetM(plhs[1], c_size);
	mxSetN(plhs[1], 1);
	return c_size;
    
    

    
    



}

//!Assumes both are sorted
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{



    double  a_subs[3]={0,0,0};
    double  b_subs[3]={0,0,0};
    mwSize a_data_length;
    mwSize b_data_length;

    

    //mexPrintf("entered\n");
	double op;
	mxClassID a_id = check_sparse_params_merge( nrhs,prhs, &a_data_length, & b_data_length, &op);
    //mexPrintf("checked params\n");
    switch (a_id)
    {
    case mxDOUBLE_CLASS:

        perform_merge(plhs,a_id,(double *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]),(double *)mxGetData(prhs[2]),(index_type *)mxGetData(prhs[3]),a_data_length,b_data_length,op);
        break;
    case mxSINGLE_CLASS:
		perform_merge(plhs, a_id, (float *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (float *)mxGetData(prhs[2]), (index_type *)mxGetData(prhs[3]), a_data_length, b_data_length, op);
        break;
    case mxINT32_CLASS:
		perform_merge(plhs, a_id, (int32_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (int32_t *)mxGetData(prhs[2]), (index_type *)mxGetData(prhs[3]), a_data_length, b_data_length, op);
        break;
	case mxINT64_CLASS:
		perform_merge(plhs, a_id, (int64_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (int64_t *)mxGetData(prhs[2]), (index_type *)mxGetData(prhs[3]), a_data_length, b_data_length, op);
		break;
	case mxUINT32_CLASS:
		perform_merge(plhs, a_id, (uint32_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (uint32_t *)mxGetData(prhs[2]), (index_type *)mxGetData(prhs[3]), a_data_length, b_data_length, op);
		break;
	case mxUINT64_CLASS:
		perform_merge(plhs, a_id, (uint64_t *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (uint64_t *)mxGetData(prhs[2]), (index_type *)mxGetData(prhs[3]), a_data_length, b_data_length, op);
		break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



}

