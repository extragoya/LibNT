#include "mex.h"

#include "DenseLattice.h"
#include "MappedDenseLattice.h"
#include "Util.h"
#include "LatticeException.h"




template<class T>
void perform_mult(T * C,T*A,  T*B, const mwSize*a_subs,const mwSize*b_subs){


    const LibMIA::MappedDenseLattice<T> latA(A,a_subs[0],a_subs[1],a_subs[2]);
    const LibMIA::MappedDenseLattice<T> latB(B,b_subs[0],b_subs[1],b_subs[2]);
    mwSize c_subs[]={a_subs[0], b_subs[1], a_subs[2]};

    LibMIA::MappedDenseLattice<T> latC(C,c_subs[0],c_subs[1],c_subs[2]);

    try{
        latC=latA*latB;
    }
    catch(LibMIA::LatticeException& e){
        mexErrMsgTxt(e.what());

    }



}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


    if (nrhs!=2)
        mexErrMsgTxt("Two arguments must be provided to dense lattice multiplication.");
    else if(nlhs>1)
        mexErrMsgTxt("Too many output arguments.");
    else if(nlhs==0)
        mexErrMsgTxt("No output");
    bool _valid= mxIsNumeric (prhs[0])& mxIsNumeric (prhs[1]);

    if (!_valid)
        mexErrMsgTxt("Invalid class\n");

    mxClassID a_id =mxGetClassID(prhs[0]);
    mxClassID b_id =mxGetClassID(prhs[1]);
    if (a_id!=b_id)
        mexErrMsgTxt("Lattice data must have the same underlying data type");



    mwSize a_subs_l=mxGetNumberOfDimensions(prhs[0]);
    mwSize b_subs_l=mxGetNumberOfDimensions(prhs[1]);
    if (a_subs_l!=b_subs_l)
        mexErrMsgTxt("Lattice data must have the same number of dimensions");
    if (a_subs_l!=3 && a_subs_l!=2)
        mexErrMsgTxt("Lattice data must be three-dimensional");

    const mwSize* a_subs=mxGetDimensions(prhs[0]);
    const mwSize* b_subs=mxGetDimensions(prhs[1]);

    mwSize c_subs[]={a_subs[0], b_subs[1], 1};
    mwSize a_subs_m[]={a_subs[0], a_subs[1],1};
    mwSize b_subs_m[]={b_subs[0], b_subs[1],1};
    if (a_subs_l==3){
        a_subs_m[2]=a_subs[2];
        b_subs_m[2]=b_subs[2];
        c_subs[2]=a_subs[2];
    }



    plhs[0]=mxCreateNumericArray(3, c_subs,   a_id, mxREAL);
    //plhs[0]=mxCreateNumericArray(0, 0,   a_id, mxREAL);

    switch (a_id)
    {
    case mxDOUBLE_CLASS:

        perform_mult((double *)mxGetData(plhs[0]),(double *)mxGetData(prhs[0]), (double *)mxGetData(prhs[1]),a_subs_m,b_subs_m);
        break;
    case mxSINGLE_CLASS:
        perform_mult((float *)mxGetData(plhs[0]),(float *)mxGetData(prhs[0]), (float *)mxGetData(prhs[1]),a_subs_m,b_subs_m);
        break;
    case mxINT32_CLASS:
        perform_mult((int *)mxGetData(plhs[0]),(int *)mxGetData(prhs[0]), (int *)mxGetData(prhs[1]),a_subs_m,b_subs_m);
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
