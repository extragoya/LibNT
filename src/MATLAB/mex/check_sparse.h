#ifndef CHECK_SPARSE_H_INCLUDED
#define CHECK_SPARSE_H_INCLUDED
#include <algorithm>
#include <type_traits>
#include <stdint.h>

#include "mex.h"
typedef int64_t mex_index_type;
typedef std::make_unsigned<mex_index_type>::type mex_unsigned_type;
mxClassID check_sparse_params_lattice(int nrhs, const mxArray *prhs[], mex_index_type  *a_subs[], mex_index_type  *b_subs[], mwSize* a_data_length, mwSize* b_data_length, bool allowed_seven=false);

mxClassID check_sparse_params_merge(int nrhs, const mxArray *prhs[], mwSize* a_data_length, mwSize* b_data_length,double *op);

mxClassID check_sparse_unary(int nrhs, const mxArray *prhs[], mwSize* a_data_length);

mwSize check_sparse_indices(int nrhs, const mxArray *prhs[],int param_index );

mwSize check_sparse_linIdx(int nrhs, const mxArray *prhs[], int param_index);

mwSize check_sparse_dims(int nrhs, const mxArray *prhs[], int param_index);

mxClassID check_sparse_data(int nrhs, const mxArray *prhs[], int param_index, mwSize* data_length);

double check_passed_double(int nrhs, const mxArray *prhs[], int param_index);

mxClassID check_dense_params_lattice(int nrhs, const mxArray *prhs[], mwSize  *a_subs, int param_index);

#endif // CHECK_SPARSE_H_INCLUDED
