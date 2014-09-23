#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdint>
#include <numeric>
#include "SparsePermuteMex_Export.h" //must include before mex.h so mex.h can use macro definitions
#include "mex.h"
#include "IndexUtil.h"
#include "LibMIARadix.h"
#include "check_sparse.h"


template<class T>
void perform_permute(T*a_data, mex_index_type  * a_indices, mwSize a_size, const std::vector<unsigned char> & curIdx, const std::vector<unsigned char> & desiredIdx, const std::vector<mex_index_type> & dims){
	using namespace LibMIA::internal;
	

	auto reverseShuffleSequence = getShuffleSequence(curIdx, desiredIdx); //get the shuffle sequence from new to old
	std::vector<mex_unsigned_type> divisors;
	std::vector<mex_unsigned_type> max_sizes;
	auto shuffled_dims = dims;
	reorder_from(dims, desiredIdx, shuffled_dims);
	bool first_stage = setupPermute(reverseShuffleSequence, shuffled_dims, divisors, max_sizes);
	auto dimensionality = std::accumulate(dims.begin(), dims.end(), 1);
	//create RadixShuffle object
	RadixShuffle<mex_index_type, T, 2048, 11, 3000> radixShuffle(max_sizes, divisors, dimensionality, first_stage);
	//permute the sparse data based on the stage information provided
	radixShuffle.permute(a_indices, a_data, a_size);

}



//!assumes indices have already been changed to desired lexicographical order
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	
	if (nrhs != 5)
		mexErrMsgTxt("Five arguments must be provided to changeLinIdx function.");
	mwSize a_data_length;
	mxClassID a_id = check_sparse_data(nrhs, prhs, 0,&a_data_length);
	mwSize a_index_length = check_sparse_indices(nrhs, prhs, 1);
	mwSize a_degree = check_sparse_linIdx(nrhs, prhs, 2);
	mwSize a_degree2 = check_sparse_dims(nrhs, prhs, 4);

	if (a_data_length != a_index_length)
		mexErrMsgTxt("Vector of data and indices must be the same length.");
	
	if (a_degree!=a_degree2)
		mexErrMsgTxt("Lexicographical indices vectors and dimensions vectors must be same size.");
	if (!a_data_length)
		return;
	
	mex_index_type* a_indices = (mex_index_type *)mxGetData(prhs[1]);
	unsigned char * current_permute_idx = (unsigned char *)mxGetData(prhs[2]);
	unsigned char * desired_permute_idx = (unsigned char *)mxGetData(prhs[3]);
	mex_index_type* _dims = (mex_index_type *)mxGetData(prhs[4]);

	std::vector<unsigned char> curIdx(a_degree);
	std::vector<unsigned char> desiredIdx(a_degree);
	std::vector<mex_index_type> dims(a_degree);
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
	
	switch (a_id)
	{
	case mxDOUBLE_CLASS:

		perform_permute((double *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), a_data_length, curIdx, desiredIdx, dims);
		break;
	case mxSINGLE_CLASS:
		perform_permute((float *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), a_data_length, curIdx, desiredIdx, dims);
		
		break;
	case mxINT32_CLASS:
		perform_permute((int32_t *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), a_data_length, curIdx, desiredIdx, dims);
		
		break;
	case mxINT64_CLASS:
		perform_permute((int64_t *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), a_data_length, curIdx, desiredIdx, dims);
		
		break;
	case mxUINT32_CLASS:
		perform_permute((uint32_t *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), a_data_length, curIdx, desiredIdx, dims);
		
		break;
	case mxUINT64_CLASS:
		perform_permute((uint64_t *)mxGetData(prhs[0]), (mex_index_type *)mxGetData(prhs[1]), a_data_length, curIdx, desiredIdx, dims);
		
		break;
	case mxUNKNOWN_CLASS:
		mexErrMsgTxt("Unrecognized data type");
		break;
	default:
		mexErrMsgTxt("Data type must be arithmetic");
		break;
	}



}

