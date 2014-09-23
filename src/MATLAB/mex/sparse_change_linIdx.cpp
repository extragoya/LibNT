#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdint>
#include "SparseChangeLinIdxMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "IndexUtil.h"
#include "check_sparse.h"







void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	using namespace LibMIA::internal;
	if (nrhs != 4)
		mexErrMsgTxt("Four arguments must be provided to changeLinIdx function.");
	
	mwSize a_data_length = check_sparse_indices(nrhs, prhs,0);
	mwSize a_degree = check_sparse_linIdx(nrhs, prhs, 1);
	mwSize a_degree2 = check_sparse_dims(nrhs, prhs, 3);

	if (a_degree!=a_degree2)
		mexErrMsgTxt("Lexicographical indices vectors and dimensions vectors must be same size.");
	if (!a_data_length)
		return;
	
	mex_index_type* a_indices = (mex_index_type *)mxGetData(prhs[0]);
	unsigned char * current_permute_idx = (unsigned char *)mxGetData(prhs[1]);
	unsigned char * desired_permute_idx = (unsigned char *)mxGetData(prhs[2]);
	mex_index_type* _dims = (mex_index_type *)mxGetData(prhs[3]);

	std::vector<unsigned char> curIdx(a_degree);
	std::vector<unsigned char> desiredIdx(a_degree);
	std::vector<unsigned char> dims(a_degree);
	std::copy(current_permute_idx, current_permute_idx + a_degree, curIdx.begin());
	std::copy(desired_permute_idx, desired_permute_idx + a_degree, desiredIdx.begin());
	std::copy(_dims, _dims + a_degree, dims.begin());
	//set indices to be zero-indexed
	for (auto &i : curIdx)
		i--;
	for (auto &i : desiredIdx)
		i--;
	if (std::equal(curIdx.begin(),curIdx.end(),desiredIdx.begin()))
		return;
	
	auto reorder_Dims = dims;
	auto new_reorder_Dims = dims;

	reorder_from(dims, curIdx, reorder_Dims);
	reorder_from(dims, desiredIdx, new_reorder_Dims);
	auto dim_accumulator = createDimAccumulator(reorder_Dims);
	auto dim_accumulator_fast = createDimAccumulator_libdivide(reorder_Dims);
	auto multiplier = createMultiplier(new_reorder_Dims);
	auto index_order = getShuffleSequence(desiredIdx, curIdx);
	auto new_multiplier = multiplier;
	reorder_from(multiplier, index_order, new_multiplier);

	for (size_t i = 0; i < a_data_length;++i){
		a_indices[i] = reShuffleLinearIndex(a_indices[i], new_multiplier, dim_accumulator_fast, dim_accumulator);
	}



}

