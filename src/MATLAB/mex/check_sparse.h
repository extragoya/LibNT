#ifndef CHECK_SPARSE_H_INCLUDED
#define CHECK_SPARSE_H_INCLUDED
#include <algorithm>
#include "mex.h"

mxClassID check_sparse_params(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[],double  *a_subs,double  *b_subs,mwSize* a_data_length,mwSize* b_data_length);

#endif // CHECK_SPARSE_H_INCLUDED
