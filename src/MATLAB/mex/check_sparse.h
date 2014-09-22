#ifndef CHECK_SPARSE_H_INCLUDED
#define CHECK_SPARSE_H_INCLUDED
#include <algorithm>
#include "mex.h"

mxClassID check_sparse_params_lattice(int nrhs, const mxArray *prhs[],double  *a_subs,double  *b_subs,mwSize* a_data_length,mwSize* b_data_length);

mxClassID check_sparse_params_merge(int nrhs, const mxArray *prhs[], mwSize* a_data_length, mwSize* b_data_length,double *op);

#endif // CHECK_SPARSE_H_INCLUDED
