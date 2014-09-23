#ifndef CHECK_SPARSE_H_INCLUDED
#define CHECK_SPARSE_H_INCLUDED
#include <algorithm>
#include <type_traits>

#include "mex.h"
typedef long long mex_index_type;
typedef std::make_unsigned<mex_index_type>::type mex_unsigned_type;
mxClassID check_sparse_params_lattice(int nrhs, const mxArray *prhs[],double  *a_subs,double  *b_subs,mwSize* a_data_length,mwSize* b_data_length);

mxClassID check_sparse_params_merge(int nrhs, const mxArray *prhs[], mwSize* a_data_length, mwSize* b_data_length,double *op);

mxClassID check_sparse_unary(int nrhs, const mxArray *prhs[], mwSize* a_data_length);

mwSize check_sparse_indices(int nrhs, const mxArray *prhs[],int param_index );

mwSize check_sparse_linIdx(int nrhs, const mxArray *prhs[], int param_index);

mwSize check_sparse_dims(int nrhs, const mxArray *prhs[], int param_index);

mxClassID check_sparse_data(int nrhs, const mxArray *prhs[], int param_index, mwSize* data_length);

#endif // CHECK_SPARSE_H_INCLUDED
