#include <iostream>
#include <sstream>
#include <string>
#include "SparseLatticeMultMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"
#include "SparseLattice.h"
#include "MappedSparseLattice.h"
#include "LibMIAException.h"




template<class T>
mwSize perform_mult(mxArray *plhs[], mxClassID a_id, T*a_data, mex_index_type  * a_indices, T*b_data, mex_index_type * b_indices, mex_index_type *a_subs, mwSize a_size, mex_index_type *b_subs, mwSize b_size,int algorithm){

    T * c_data;
	mex_index_type * c_indices;

    

    try{
		LibMIA::SparseLattice<T> latC; 
		switch (algorithm){
		case 0:{
			LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, a_subs[0], a_subs[1], a_subs[2]);
			LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, b_subs[0], b_subs[1], b_subs[2]);
			latC = latA*latB;
			}
			break;
		case 1: //CSC
		{
			LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, a_subs[0], a_subs[1], a_subs[2]);
			LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, b_subs[0], b_subs[1], b_subs[2]);
			latC = latA.template csc_times<false>(latB);
		}
			break;
		case 2: //DCSC
		{
			LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, a_subs[0], a_subs[1], a_subs[2]);
			LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, b_subs[0], b_subs[1], b_subs[2]);
			latC = latA.template csc_times<true>(latB);
		}
			break;
		case 3: //CSCNA
		{
			LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, a_subs[0], a_subs[1], a_subs[2]);
			LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, b_subs[0], b_subs[1], b_subs[2]);
			latC = latA.template csc_no_accum<false>(latB);
		}
			break;
		case 4: //outer-product
		{
			LibMIA::MappedSparseLattice<T> latA(a_data, a_indices, a_size, a_subs[0], a_subs[1], a_subs[2]);
			LibMIA::MappedSparseLattice<T> latB(b_data, b_indices, b_size, b_subs[0], b_subs[1], b_subs[2], LibMIA::RowMajor); //your b data must be in row major
			latC = latA.outer_times(latB);
		}
			break;
		default:
			mexErrMsgTxt("Must specify a multiplication algorithm between 0 and 4.");
		}
        c_data = (T *)mxCalloc(latC.size() , sizeof(T));
		c_indices = (mex_index_type *)mxCalloc(latC.size(), sizeof(mex_index_type));
        std::copy(latC.data_begin(),latC.data_end(),&c_data[0]);
        std::copy(latC.index_begin(),latC.index_end(),&c_indices[0]);
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


	mex_index_type *a_subs;
	mex_index_type *b_subs;
    
    mwSize a_data_length;
    mwSize b_data_length;
	if (nlhs != 2)
		mexErrMsgTxt("Two output arguments, c_data c_indices, required to multiply SparseLattices.");
	else if (nlhs == 0)
		mexErrMsgTxt("No output");
    mxClassID a_id=check_sparse_params_lattice(nrhs, prhs,&a_subs,&b_subs,&a_data_length,&b_data_length,true);
	int algorithm = 0;
	if (nrhs == 7){
		algorithm=(int)check_passed_double(nrhs, prhs, 6);
	}


    switch (a_id)
    {
    case mxDOUBLE_CLASS:

		perform_mult(plhs, a_id, (double *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), (double *)mxGetData(prhs[3]), (mex_index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length,algorithm);
        break;
    case mxSINGLE_CLASS:
		perform_mult(plhs, a_id, (float *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), (float *)mxGetData(prhs[3]), (mex_index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length, algorithm);
        break;
    case mxINT32_CLASS:
		perform_mult(plhs, a_id, (int *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), (int *)mxGetData(prhs[3]), (mex_index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length, algorithm);
        break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized lattice data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



}

