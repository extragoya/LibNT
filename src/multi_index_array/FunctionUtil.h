// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef FUNCTION_UTIL_H
#define FUNCTION_UTIL_H
#define LIBMIA_LOG2E 1.44269504088896340736
#include <chrono>
#include <iostream>
#include <algorithm>
#include <boost/numeric/conversion/converter.hpp>

#include "LibMIAUtil.h"
#include "LibMIARanges.h"
#include "IndexUtil.h"
#include "PermuteIterator.h"
namespace LibMIA
{
#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))
namespace internal{
/** \addtogroup util Utilities
 *  @{
*/

//index_order is order going from LHS to in RHS - ie {1 3 2 0} means a(i,j,k,l)=b(l,i,k,j)
template<class Derived,class otherDerived, class Op, class index_param_type>
typename MIAMergeReturnType<Derived,otherDerived>::type
perform_implicit_merge(const MIA<Derived>& a, const MIA<otherDerived>& b,const Op& op,const std::array<index_param_type,internal::order<Derived>::value>& index_order){


    typedef typename MIAMergeReturnType<Derived,otherDerived>::type retMIAType;



    typedef typename internal::index_type<otherDerived>::type b_index_type;
    retMIAType c(a.dims());




    typename MIA<Derived>::accumulator_type dim_accumulator;
    typename MIA<Derived>::fast_accumulator_type fast_dim_accumulator;
    typename MIA<Derived>::multiplier_type multiplier;
    internal::create_shuffle_needs(a.dims(),b.dims(),index_order,dim_accumulator,fast_dim_accumulator,multiplier);
    if(a.dimensionality()>=PARALLEL_TOL){
        #pragma omp parallel for
        for(b_index_type idx=0;idx<a.dimensionality();++idx)
            c.atIdx(idx)=op(a.atIdx(idx),a.convert(b.atIdx(internal::reShuffleLinearIndex(idx,multiplier,fast_dim_accumulator,dim_accumulator))));
    }
    else{

        for(b_index_type idx=0;idx<a.dimensionality();++idx)
            c.atIdx(idx)=op(a.atIdx(idx),a.convert(b.atIdx(internal::reShuffleLinearIndex(idx,multiplier,fast_dim_accumulator,dim_accumulator))));
    }
    return c;

}


template<class Derived,class otherDerived, class Op>
typename MIAMergeReturnType<Derived,otherDerived>::type
perform_implicit_merge(const MIA<Derived>& a, const MIA<otherDerived>& b,const Op& op){



    typedef typename MIAMergeReturnType<Derived,otherDerived>::type retMIAType;
    retMIAType c(a.dims());

    if(a.dimensionality()>=PARALLEL_TOL){
        #pragma omp parallel for
        for(size_t idx=0;idx<a.dimensionality();++idx)
            c.atIdx(idx)=op(a.atIdx(idx),a.convert(b.atIdx(idx)));
    }
    else{
        #pragma omp parallel for
        for(size_t idx=0;idx<a.dimensionality();++idx)
            c.atIdx(idx)=op(a.atIdx(idx),a.convert(b.atIdx(idx)));
    }
    return c;

//    typedef typename internal::function_type<retMIAType>::type function_type;
//    typedef typename internal::index_type<retMIAType>::type index_type;
//
//
//
//
//    auto _function=[&a,&b,op](index_type idx){
//            return op(a.atIdx(idx),a.convert(b.atIdx(idx)));
//            //return a.atIdx(idx)+b.atIdx(idx);
//    };
//
//
//    return retMIAType(_function,a.dims());

}






//Base Case
template<size_t cur_partition,class MIAType2,class index_it,size_t no_con_partition,
    typename boost::enable_if_c< cur_partition == 0, int >::type = 0
    >
typename internal::data_type<MIAType2>::type collect_contract_partitions(const MIAType2 & source,typename internal::index_type<MIAType2>::type cur_index,const std::array<size_t, no_con_partition>& contract_ranges,const index_it contract_idx_end,const std::array<int,no_con_partition> & contract_partition)
{


    return source.atIdx(cur_index);


}

template<size_t cur_partition,class MIAType2,class index_it,size_t no_con_partition,
    typename boost::disable_if_c< cur_partition == 0, int >::type = 0
    >
typename internal::data_type<MIAType2>::type collect_contract_partitions(const MIAType2 & source,typename internal::index_type<MIAType2>::type cur_index,const std::array<size_t, no_con_partition>& contract_ranges,const index_it contract_idx_end,const std::array<int,no_con_partition> & contract_partition)
{

    typedef typename internal::data_type<MIAType2>::type data_type;
    data_type sum=0;
    for(int j=0;j<(int)contract_ranges[cur_partition-1];++j){
        auto j_contract_idx=internal::get_contract_idx(j, contract_idx_end-contract_partition[cur_partition-1],contract_idx_end, source.dims());
        sum+=collect_contract_partitions<cur_partition-1,MIAType2,index_it,no_con_partition>(source,cur_index+j_contract_idx,contract_ranges,contract_idx_end-contract_partition[cur_partition-1],contract_partition);

    }
    return sum;


}

template<class Operand1, class Operand2, class index_type, size_t degree>
index_type performShuffledOperation(Operand1 & restrict_libmia operand1, const Operand2 & restrict_libmia operand2, const std::array<index_type, degree> & dims, const std::array<index_type, degree> & multiplier, size_t curIndex, index_type sourceIdx, index_type destIdx){
	if (curIndex == 0){
		const index_type end = destIdx + dims[0];
		
		
			
		for (; destIdx < end; ++destIdx,sourceIdx+=multiplier[0]){
			operand1.atIdx(destIdx) = operand2.atIdx(sourceIdx);
				
		}
		
	}		
	else{
		
		for (index_type dim_idx = 0; dim_idx < dims[curIndex]; ++dim_idx){
			destIdx = performShuffledOperation(operand1, operand2, dims, multiplier, curIndex - 1, sourceIdx, destIdx);
			sourceIdx += multiplier[curIndex];
		}
	}
	return destIdx;
}

template<class Operand1, class Operand2, class index_type, size_t degree>
index_type performShuffledOperationParallel(Operand1 & restrict_libmia operand1, const Operand2 & restrict_libmia operand2, const std::array<index_type, degree> & dims, const std::array<index_type, degree> & multiplier, size_t curIndex, index_type sourceIdx, index_type destIdx){
	
	if (curIndex == 0){
		const index_type end = destIdx + dims[0];
		
		
		#pragma omp parallel for
		for (index_type temp_idx = 0; temp_idx<dims[0]; ++temp_idx){
			operand1.atIdx(destIdx + temp_idx) = operand2.atIdx(sourceIdx + temp_idx*multiplier[0]);
			
		}
		
	}
	else{
		const index_type multiplier_dest = operand1.dimensionality() / dims.back();
		#pragma omp parallel for
		for (index_type dim_idx = 0; dim_idx < dims[curIndex]; ++dim_idx){
			performShuffledOperation(operand1, operand2, dims, multiplier, curIndex - 1, sourceIdx + dim_idx*multiplier[curIndex], destIdx + multiplier_dest*dim_idx);
			//sourceIdx += multiplier[curIndex];
		}
	}
	return destIdx;
}

template<typename Derived, class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
auto latticeCopy(const MIA<Derived> &mia, const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices)
->DenseLattice<typename internal::data_type<Derived>::type>
{

	using namespace std::chrono;
	//typedef std::chrono::duration<float> float_seconds;
	high_resolution_clock::time_point t1, t2;
	//t1 = high_resolution_clock::now();
    //print_array(row_indices,"row_indices");
    //print_array(column_indices,"column_indices");
    //print_array(tab_indices,"tab_indices");
	typedef typename MIA<Derived>::unsigned_index_type unsigned_index_type;
    typedef typename internal::data_type<Derived>::type data_type;
    constexpr auto order= internal::order<Derived>::value ;
	static_assert(internal::check_index_compatibility<unsigned_index_type, idx_typeR>::type::value, "Must use an array convertable to index_type");
	static_assert(internal::check_index_compatibility<unsigned_index_type, idx_typeC>::type::value, "Must use an array convertable to index_type");
	static_assert(internal::check_index_compatibility<unsigned_index_type, idx_typeT>::type::value, "Must use an array convertable to index_type");
	static_assert(R + C + T == order, "Size of all three arrays must equal mOrder");
    //statically check number of indices match up
    size_t row_size=1, column_size=1, tab_size=1;

    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    row_size=internal::dimensionality_from(mia.dims(), row_indices);
    column_size=internal::dimensionality_from(mia.dims(), column_indices);
    tab_size=internal::dimensionality_from(mia.dims(), tab_indices);

	std::array<unsigned_index_type, R + C + T> shuffled_dims;
    std::array<size_t,R+C+T> index_order;
    concat_arrays(row_indices,column_indices,tab_indices,index_order);
    internal::reorder_from(mia.dims(),index_order,shuffled_dims);

    //std::cout<< "Tab dims " << tab_dims[0] << "\n";
    DenseLattice<data_type> lat(row_size, column_size, tab_size);



    typename MIA<Derived>::accumulator_type dim_accumulator;
    typename MIA<Derived>::fast_accumulator_type fast_dim_accumulator;
    typename MIA<Derived>::multiplier_type multiplier;
    internal::create_shuffle_needs(shuffled_dims,mia.dims(),index_order,dim_accumulator,fast_dim_accumulator,multiplier);

	if (mia.dimensionality() >= PARALLEL_TOL)
		performShuffledOperationParallel(lat, mia, shuffled_dims, multiplier, size_t(R + C + T - 1), unsigned_index_type(0), unsigned_index_type(0));
	else
		performShuffledOperation(lat, mia, shuffled_dims, multiplier, size_t(R + C + T - 1), unsigned_index_type(0), unsigned_index_type(0));

    /*if(mia.dimensionality()>=PARALLEL_TOL){
        #pragma omp parallel for
        for(index_type idx=0;idx<mia.dimensionality();++idx){
            lat.atIdx(idx)=mia.atIdx(internal::reShuffleLinearIndex(idx,multiplier,fast_dim_accumulator,dim_accumulator));
        }
    }
    else{
        for(index_type idx=0;idx<mia.dimensionality();++idx){
            lat.atIdx(idx)=mia.atIdx(internal::reShuffleLinearIndex(idx,multiplier,fast_dim_accumulator,dim_accumulator));
        }

    }*/


	/*t2 = high_resolution_clock::now();
	std::cout << "\t" << duration_cast<float_seconds>(t2 - t1).count();*/


    return lat;


}


template<typename MIA,typename otherMIA, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto implicitNoLatticeMult(const MIA &a,const otherMIA &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<MIA,otherMIA,L_outer+R_outer+Inter>::type
{

    typedef typename MIANoLatticeProductReturnType<MIA,otherMIA,L_outer+R_outer+Inter>::type RetType;
    typedef typename internal::index_type<RetType>::type index_type;
    typedef typename internal::index_type<MIA>::type a_index_type;
    typedef typename internal::index_type<otherMIA>::type b_index_type;
    typedef typename internal::function_type<RetType>::type function_type;
    static_assert(internal::check_index_compatibility<index_type,array_type>::type::value,"Must use an array convertable to index_type");
    std::array<a_index_type, Inter> l_inter_dims;
    std::array<a_index_type, L_outer> l_outer_dims;
    std::array<b_index_type, Inter> r_inter_dims;
    std::array<b_index_type, R_outer> r_outer_dims;




    //get inter and outer dimensionality and the individual dimensions that make up that number - should default to one if any of the arrays are empty
    size_t l_inter_size=internal::reorder_from(a.dims(), l_inter_idx,l_inter_dims);
    size_t r_inter_size= internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
    if(l_inter_size!=r_inter_size || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");
    internal::reorder_from(a.dims(), l_outer_idx,l_outer_dims);
    internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);

    std::array<index_type,L_outer+R_outer+Inter> retDims;
    concat_arrays(l_outer_dims, r_outer_dims,l_inter_dims,retDims);
    RetType c(retDims);
    //create lambda function
    function_type _function=[&a,&b,&c,l_inter_idx,r_inter_idx,l_outer_idx,r_outer_idx](index_type _index){
        auto full_indices=internal::ind2sub(_index,c.dims());
        a_index_type l_idx(0);
        b_index_type r_idx(0);
        l_idx+=internal::sub2ind(full_indices.begin(),full_indices.begin()+L_outer,l_outer_idx,a.dims());
        l_idx+=internal::sub2ind(full_indices.begin()+L_outer+R_outer,full_indices.end(),l_inter_idx,a.dims());
        r_idx+=internal::sub2ind(full_indices.begin()+L_outer,full_indices.begin()+L_outer+R_outer,r_outer_idx,b.dims());
        r_idx+=internal::sub2ind(full_indices.begin()+L_outer+R_outer,full_indices.end(),r_inter_idx,b.dims());
//        print_array(full_indices,"full_indices");
//        print_array(l_outer_idx,"l_outer_idx");
//        std::cout << "Outer l_idx " << internal::sub2ind(full_indices.begin(),full_indices.begin()+L_outer,l_outer_idx,a.dims()) << std::endl;
//        std::cout << "L_outer " << L_outer << " R_outer " << R_outer << " Inter " << Inter << std::endl;
//        std::cout << "l_idx " << l_idx << " r_idx " << r_idx << std::endl;
        return a.atIdx(l_idx)*b.atIdx(r_idx);
    };
    c.get_function()=_function;
    return c;



}



template<typename MIA,typename otherMIA, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto noLatticeMult(const MIA &a,const otherMIA &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<MIA,otherMIA,L_outer+R_outer+Inter>::type
{

    typedef typename MIANoLatticeProductReturnType<MIA,otherMIA,L_outer+R_outer+Inter>::type RetType;
    typedef typename internal::index_type<RetType>::type index_type;
    typedef typename internal::index_type<MIA>::type a_index_type;
    typedef typename internal::index_type<otherMIA>::type b_index_type;

    static_assert(internal::check_index_compatibility<index_type,array_type>::type::value,"Must use an array convertable to index_type");
    std::array<a_index_type, Inter> l_inter_dims;
    std::array<a_index_type, L_outer> l_outer_dims;
    std::array<b_index_type, Inter> r_inter_dims;
    std::array<b_index_type, R_outer> r_outer_dims;




    //get inter and outer dimensionality and the individual dimensions that make up that number - should default to one if any of the arrays are empty
    size_t l_inter_size=internal::reorder_from(a.dims(), l_inter_idx,l_inter_dims);
    size_t r_inter_size= internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
    if(l_inter_size!=r_inter_size || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");
    internal::reorder_from(a.dims(), l_outer_idx,l_outer_dims);
    internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);

    std::array<index_type,L_outer+R_outer+Inter> retDims;
    concat_arrays(l_outer_dims, r_outer_dims,l_inter_dims,retDims);
    RetType c(retDims);

    for(index_type idx=0;idx<c.dimensionality();++idx){
        auto full_indices=c.ind2sub(idx);
        a_index_type l_idx(0);
        b_index_type r_idx(0);
        l_idx+=internal::sub2ind(full_indices.begin(),full_indices.begin()+L_outer,l_outer_idx,a.dims());
        l_idx+=internal::sub2ind(full_indices.begin()+L_outer+R_outer,full_indices.end(),l_inter_idx,a.dims());
        r_idx+=internal::sub2ind(full_indices.begin()+L_outer,full_indices.begin()+L_outer+R_outer,r_outer_idx,b.dims());
        r_idx+=internal::sub2ind(full_indices.begin()+L_outer+R_outer,full_indices.end(),r_inter_idx,b.dims());
        c.atIdx(idx)=a.atIdx(l_idx)*b.atIdx(r_idx);


    }


    return c;



}

//assumes C, A, and B are of the same dimensions and in the same sort order (and A and B are sorted), should work for either SparseLattices or SparseMIAs
template<class C_Class, class B_Class, class A_Class,class Op>
void outside_merge_sparse_storage_containers(C_Class & C, const A_Class & A, const B_Class & B, Op op)
{

	using namespace boost::numeric;
    typedef typename internal::data_type<A_Class>::type a_data_type;

    C.clear();
    C.reserve(A.size()+B.size());
    auto a_begin=A.index_begin();
    auto b_begin=B.index_begin();
    auto a_end=A.index_end();
    auto b_end=B.index_end();
    while(a_begin<a_end && b_begin<b_end){
        if (*a_begin<*b_begin){
            C.push_back(C.convert(A.data_at(a_begin)),*a_begin);
            a_begin++;
        }
        else if  (*b_begin<*a_begin){
            C.push_back(C.convert(op(a_data_type(0),B.data_at(b_begin))),*b_begin);
            b_begin++;
        }
        else{
            C.push_back(C.convert(op(A.data_at(a_begin),B.data_at(b_begin))),*a_begin);
            a_begin++;
            b_begin++;
        }

    }
    if (a_begin==a_end){
        while (b_begin<b_end){
            C.push_back(C.convert(op(a_data_type(0),B.data_at(b_begin))),*b_begin);
            b_begin++;
        }
    }
    else{
        while (a_begin<a_end){
            C.push_back(C.convert(A.data_at(a_begin)),*a_begin);
            a_begin++;
        }
    }




}



////must be boost::tuples of iterators. Assumes a's container is sized to be a.size+b.size
//template<class AStorageItType, class BStorageItType, class Op>
//AStorageItType merge_sparse_storage_containers(AStorageItType  a_begin,AStorageItType  a_end,BStorageItType  b_begin,BStorageItType  b_end,Op op)
//{
//    using namespace boost::numeric;
//    typedef typename boost::remove_reference<typename BStorageItType::value_type::first_type>::type b_data_type;
//    typedef typename boost::remove_reference<typename AStorageItType::value_type::first_type>::type a_data_type;
//
//    typedef converter<a_data_type,b_data_type,conversion_traits<a_data_type,b_data_type>,def_overflow_handler,RoundEven<b_data_type>> to_mdata_type;
//    AStorageItType a_actual_end=a_end;
//    AStorageItType a_actual_begin=a_begin;
//    while(a_begin<a_end && b_begin<b_end){
//        if (std::get<1>(*a_begin)<std::get<1>(*b_begin)){
//            a_begin++;
//        }
//        else if  (std::get<1>(*b_begin)<std::get<1>(*a_begin)){
//            std::get<0>(*a_actual_end)=op(a_data_type(0),to_mdata_type::convert(std::get<0>(*b_begin)));
//            std::get<1>(*a_actual_end++)=std::get<1>(*b_begin++);
//
//
//        }
//        else{
//            std::get<0>(*a_begin)=op(std::get<0>(*a_begin),to_mdata_type::convert(std::get<0>(*b_begin)));
//            a_begin++;
//            b_begin++;
//        }
//
//    }
//    if (a_begin==a_end){
//        while (b_begin<b_end){
//            std::get<0>(*a_actual_end)=op(a_data_type(0),to_mdata_type::convert(std::get<0>(*b_begin)));
//            std::get<1>(*a_actual_end++)=std::get<1>(*b_begin++);
//        }
//    }
//
//    std::inplace_merge(a_actual_begin,a_end,a_actual_end,[](const typename AStorageItType::value_type& lhs, const typename AStorageItType::value_type& rhs)
//    {
//        return std::get<1>(lhs)<std::get<1>(rhs);
//    });
//
//
//    return a_actual_end;
//
//}
//
////must be boost::tuples of iterators. Assumes a's container is sized to be a.size+b.size
//template<class ADataIt, class AIndexIt, class BDataIt, class BIndexIt,class Op>
//ADataIt merge_sparse_storage_containers(ADataIt  a_data_begin,ADataIt  a_data_end,AIndexIt  a_index_begin,AIndexIt  a_index_end,BDataIt  b_data_begin,BDataIt  b_data_end,BIndexIt  b_index_begin,BIndexIt  b_index_end,Op op)
//{
//    using namespace boost::numeric;
//    typedef typename ADataIt::value_type a_data_type;
//    typedef typename BDataIt::value_type b_data_type;
//
//
//    typedef converter<a_data_type,b_data_type,conversion_traits<a_data_type,b_data_type>,def_overflow_handler,RoundEven<b_data_type>> to_mdata_type;
//    ADataIt a_actual_data_end=a_data_end;
//    AIndexIt a_actual_index_end=a_index_end;
//    ADataIt a_cur_data_it=a_data_begin;
//    AIndexIt a_cur_index_it=a_index_begin;
//    while(a_cur_data_it<a_data_end && b_data_begin<b_data_end){
//        if (*a_cur_index_it<*b_index_begin){
//            a_cur_index_it++;
//            a_cur_data_it++;
//        }
//        else if  (*b_index_begin<*a_cur_index_it){
//            *a_actual_data_end++=*b_data_begin++;
//            *a_actual_index_end++=*b_index_begin++;
//
//
//        }
//        else{
//            *a_cur_data_it=op(*a_cur_data_it,to_mdata_type::convert(*b_data_begin++));
//            a_cur_data_it++;
//            a_cur_index_it++;
//            b_index_begin++;
//        }
//
//    }
//    if (a_cur_data_it==a_data_end){
//        while (b_data_begin<b_data_end){
//            *a_actual_data_end++=*b_data_begin++;
//            *a_actual_index_end++=*b_index_begin++;
//
//        }
//    }
////    std::cout << "Index\t Data in scan merge" << std::endl;
////    auto j=a_index_begin;
////    for(auto i=a_data_begin;i<a_actual_data_end;++i,++j)
////        std::cout << *j << "\t " << *i << std::endl;
////
////    std::cout << std::endl;
////
////    std::cout << " diff " << a_index_end-a_index_begin << std::endl;
//    std::inplace_merge(make_sort_permute_iter(a_index_begin,a_data_begin),
//                    make_sort_permute_iter(a_index_end,a_data_end),
//              make_sort_permute_iter(a_actual_index_end,a_actual_data_end),
//                    sort_permute_iter_compare<AIndexIt,ADataIt>());
//
////    std::inplace_merge(a_data_begin,a_data_end,a_actual_data_end,[&](const typename ADataIt::value_type& lhs, const typename ADataIt::value_type& rhs)
////    {
////        return *(a_index_begin+(&lhs-&(*a_data_begin))) <*(a_index_begin+(&rhs-&(*a_data_begin)));
////    });
////    std::inplace_merge(a_index_begin,a_index_end,a_actual_index_end);
//    /*std::cout << " diff " << a_index_end-a_index_begin << std::endl;
//    std::cout << "Index\t Data in AFTER scan merge" << std::endl;
//    for(auto i=a_data_begin,j=a_index_begin;i<a_actual_data_end;++i,++j)
//        std::cout << *j << "\t " << *i << std::endl;
//
//    std::cout << std::endl;*/
//    return a_actual_data_end;
//
//}


template<typename index_type,typename T,
    typename boost::enable_if<
        boost::is_integral<T>,
        int
    >::type=0
>
Range<index_type> create_range(T t){
    return Range<index_type>((index_type)t); //create range object based on t;
}

//if the current variadic template parameter is a Range object, simply just return it
template<typename index_type>
Range<index_type> & create_range(Range<index_type>&t){
    return t;
}

//base case for variadic parameters - do not alter the iterator
template<typename ItType>
void get_range_array(ItType it){
    return;
}



//if the current variable in the varadic template is an integral type, create a Range object for that type and add it using the current array iterator
template<typename ItType,typename T,typename...Ranges>
void get_range_array(ItType it,T t,Ranges...ranges){
    typedef typename ItType::value_type RangeType;
    typedef typename RangeType::index_type index_type;
    *it=create_range<index_type>(t);

    get_range_array(it+1,ranges...); //recurse to next spot in the container and the next variadic template parameter

}




//! Converts a scalar value to data_type
/*!
    \tparam from_data_type the data_type you are converting from
*/
template<class data_type,class from_data_type,typename boost::enable_if< boost::is_pod< from_data_type >, int >::type = 0>
inline data_type convert(const from_data_type from){
    using namespace boost::numeric;
    typedef boost::numeric::converter<data_type,boost::uniform_real<>::result_type> to_mdata_type;
    return to_mdata_type::convert(from);
}


struct print_class_name {
    template <typename T>
    void operator()( T t ) const {
       std::cout << typeid(t).name() << " ";
    }
};

inline long double log2(const long double x){
	return  std::log(x) * LIBMIA_LOG2E;
}



template<typename T1>
inline T1 manual_int_power(const T1 base,const int _exp){
    T1 result=1;
    for(int i=0;i<_exp;++i)
    {
        result*=base;

    }
    return result;

}

/*! @} */
}



template<class T>
struct MIAprint
{
    void operator() (T i)
    {
        std::cout << " " << i;
    }
} ;

template<class T>
struct select_first
{
    T& operator()(T&left, T& right){
        return left;
    }
};



template<class array_type>
void print_array(const array_type & _array, const std::string &header){
    std::cout << header;
    for(auto & _i:_array){
        std::cout << " " << _i;
    }
    std::cout << std::endl;

}

template<class array_type>
void print_array_on_line(const array_type & _array){
    for(auto & _i:_array){
        std::cout << " " << _i;
    }


}

template<class T1, class T2,size_t _size>
bool compare_arrays(const std::array<T1,_size> & array1, const std::array<T2,_size> & array2){
    typedef boost::numeric::converter<T1,T2> to_mdata_type;
    for(size_t i=0;i<_size;++i)
        if (array1[i]!=to_mdata_type::convert(array2[i]))
            return false;

    return true;


}

template<class data_type>
struct array_converter
{

    template<class other_data_type,size_t _size>
    static std::array<data_type,_size> convert(const std::array<other_data_type,_size> & _from)
    {
        typedef boost::numeric::converter<data_type,other_data_type> to_mdata_type;
        std::array<data_type,_size> ret;
        for(size_t i=0;i<_size;++i)
            ret[i]=to_mdata_type::convert(_from[i]);

        return ret;
    }

    template<size_t _size>
    static std::array<data_type,_size> convert(std::array<data_type,_size> & _from){
        return _from;
    }
};

//!prec must be positive
template<typename T, typename T2,typename T3>
inline bool isEqualFuzzy(T a, T2 b, T3 prec = Tolerance<T>::tolerance)
{
  if(std::abs(a) < 1 || std::abs(b) < 1)
    return std::abs(a-b)<=prec;
  else{
	  return std::abs(a - b)<= std::min(std::abs(a), std::abs(b))*prec;
	  
  }
}


//! Removes data with duplicated indices - conflicts are solved by using the collector class, ie std::plus<data_type>
template<class index_it, class data_it,class Collector>
size_t collect_duplicates_function(index_it index_begin, index_it index_end, data_it data_begin,Collector collector)
{

	if (index_begin == index_end)
		return 0;
    auto result_idx = index_begin;
    auto result_data= data_begin;
    auto first=result_idx;
    auto first_data=data_begin;

    while (++first < index_end)
    {
        ++first_data;
        if (*result_idx != *first){
            *(++result_idx)=*first;
            *(++result_data)=*first_data;
        }
        else{

            *result_data=collector(*result_data,*first_data);
        }
    }

    return result_idx-index_begin+1;

}

}


#endif // FUNCTION_UTIL_H

