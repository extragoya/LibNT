#include <iostream>
#include <sstream>
#include <string>
#include "SparseLatticeSolveMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "check_sparse.h"
#include "SparseLattice.h"
#include "MappedSparseLattice.h"
#include "LibMIAException.h"


typedef long long index_type;

template<class T>
mwSize perform_solve(mxArray *plhs[],mxClassID a_id, T*a_data, index_type  * a_indices,  T*b_data, index_type * b_indices, double *a_subs,mwSize a_size,double *b_subs,mwSize b_size){

    T * c_data;
    index_type  * c_indices;

    LibMIA::MappedSparseLattice<T> latA(a_data,a_indices,a_size,(int)a_subs[0],(int)a_subs[1],(int)a_subs[2]);
    LibMIA::MappedSparseLattice<T> latB(b_data,b_indices,b_size,(int)b_subs[0],(int)b_subs[1],(int)b_subs[2]);

    try{
        LibMIA::SparseLattice<T> latC=latA.solve(latB);

        c_data = (T *)mxCalloc(latC.size() , sizeof(T));
        c_indices = (index_type *)mxCalloc(latC.size() , sizeof(index_type));
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
		return 0; //superfluous since mexErrMsgTxt exits program, but done to avoid warnings

    }



}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


    double  a_subs[3]={0,0,0};
    double  b_subs[3]={0,0,0};
    mwSize a_data_length;
    mwSize b_data_length;
    mxClassID a_id=check_sparse_params(nlhs, plhs, nrhs, prhs,a_subs,b_subs,&a_data_length,&b_data_length);

    switch (a_id)
    {
    case mxDOUBLE_CLASS:

        perform_solve(plhs,a_id,(double *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]),(double *)mxGetData(prhs[3]),(index_type *)mxGetData(prhs[4]),a_subs,a_data_length,b_subs,b_data_length);
        break;
    case mxSINGLE_CLASS:
        perform_solve(plhs,a_id,(float *)mxGetData(prhs[0]), (index_type *)mxGetData(prhs[1]),(float *)mxGetData(prhs[3]),(index_type *)mxGetData(prhs[4]),a_subs,a_data_length,b_subs,b_data_length);
        break;
    case mxUNKNOWN_CLASS:
        mexErrMsgTxt("Unrecognized lattice data type");
        break;
    default:
        mexErrMsgTxt("Data type must be arithmetic");
        break;
    }



}

