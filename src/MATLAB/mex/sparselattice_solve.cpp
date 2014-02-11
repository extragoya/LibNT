#include <iostream>
#include <sstream>
#include <string>
#include "SparseLatticeSolveMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"
#include "DenseLattice.h"
#include "MappedSparseLattice.h"
#include "LibMIAException.h"


typedef long long index_type;

template<class T>
void perform_solve(T * C, T*a_data, index_type  * a_indices, T*b_data, index_type * b_indices, double *a_subs, mwSize a_size, double *b_subs, mwSize b_size){

    
  

    LibMIA::MappedSparseLattice<T> latA(a_data,a_indices,a_size,(int)a_subs[0],(int)a_subs[1],(int)a_subs[2]);
    LibMIA::MappedSparseLattice<T> latB(b_data,b_indices,b_size,(int)b_subs[0],(int)b_subs[1],(int)b_subs[2]);

    try{
        LibMIA::DenseLattice<T> latC=latA.solve(latB);
		std::copy(latC.data_begin(), latC.data_end(), C);
        

    }
    catch(LibMIA::LatticeException& e){
        mexErrMsgTxt(e.what());
		

    }



}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


    double  a_subs[3]={0,0,0};
    double  b_subs[3]={0,0,0};
    mwSize a_data_length;
    mwSize b_data_length;
	if (nlhs != 1)
		mexErrMsgTxt("One output argument, a dense lattice, required to solve SparseLattices.");
	
    mxClassID a_id=check_sparse_params(nrhs, prhs,a_subs,b_subs,&a_data_length,&b_data_length);
	mwSize c_subs[3];
	c_subs[0] = a_subs[1]; c_subs[1] = b_subs[1]; c_subs[2] = a_subs[2];
	
	plhs[0] = mxCreateNumericArray(3, c_subs, a_id, mxREAL);
    switch (a_id)
    {
    case mxDOUBLE_CLASS:

		perform_solve((double *)mxGetData(plhs[0]), (double *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (double *)mxGetData(prhs[3]), (index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length);
        break;
    case mxSINGLE_CLASS:
		perform_solve((float *)mxGetData(plhs[0]), (float *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]), (float *)mxGetData(prhs[3]), (index_type *)mxGetData(prhs[4]), a_subs, a_data_length, b_subs, b_data_length);
        break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized lattice data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



}

