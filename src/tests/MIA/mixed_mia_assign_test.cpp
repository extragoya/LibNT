#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE MixedMIAAssignTests
#include "MIAConfig.h"
#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "SparseMIA.h"
#include "DenseMIA.h"
#include "Index.h"





template<class _data_type>
void mult_work(size_t dim1, size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;


    LibMIA::DenseMIA<_data_type,4> original_dense_a(dim1,dim2,dim1,dim2);
    original_dense_a.randu(0,10);
    for(auto it=original_dense_a.data_begin();it<original_dense_a.data_end();++it){
        if(*it<7)
            *it=0;
    }
    LibMIA::DenseMIA<_data_type,4> dense_a=original_dense_a;
    LibMIA::DenseMIA<_data_type,4> dense_temp;
    LibMIA::SparseMIA<_data_type,4> original_sparse_a(dim1,dim2,dim1,dim2);
    LibMIA::SparseMIA<_data_type,4> sparse_a;
    LibMIA::SparseMIA<_data_type,4> sparse_temp;
    original_sparse_a.resize(original_sparse_a.dimensionality()/2);

    original_sparse_a.randu(-50,50);
    original_sparse_a.rand_indices();
    original_sparse_a.collect_duplicates();


    sparse_a=original_dense_a;

    BOOST_CHECK_MESSAGE(sparse_a==original_dense_a,std::string("Straight dense to sparse assign for ")+typeid(_data_type).name());


    sparse_a=original_sparse_a;
    dense_a=sparse_a;
    //sparse_a.print();
    //dense_a.print();
    BOOST_CHECK_MESSAGE(sparse_a==dense_a,std::string("Straight sparse to dense assign for ")+typeid(_data_type).name());


    //now try with indices
    sparse_a(i,j,k,l)=original_dense_a(i,j,k,l);

    BOOST_CHECK_MESSAGE(sparse_a==original_dense_a,std::string("Straight dense to sparse assign with indices for ")+typeid(_data_type).name());


    sparse_a=original_sparse_a;
    dense_a(i,j,k,l)=sparse_a(i,j,k,l);

    BOOST_CHECK_MESSAGE(sparse_a==dense_a,std::string("Straight sparse to dense assign with indices for ")+typeid(_data_type).name());

    sparse_temp(i,j,k,l)=original_dense_a(k,j,l,i);
    sparse_a(k,j,l,i)=sparse_temp(i,j,k,l);


    BOOST_CHECK_MESSAGE(sparse_a==original_dense_a,std::string("Dense to sparse assign with shuffled indices for ")+typeid(_data_type).name());

    dense_temp(i,j,k,l)=original_sparse_a(k,j,l,i);
    dense_a(k,j,l,i)=dense_temp(i,j,k,l);
//    original_sparse_a.print();
//    dense_a.print();

    BOOST_CHECK_MESSAGE(dense_a==original_sparse_a,std::string("Sparse to dense assign with shuffled indices for ")+typeid(_data_type).name());


}

BOOST_AUTO_TEST_CASE( MixedMIAAssignTests )
{

    //mult_work<double>(3,3);
    //mult_work<float>(3,3);
    //mult_work<int>(3,3);
    //mult_work<long>(3,3);
//
//


    mult_work<double>(8,5);
    mult_work<float>(8,5);
    mult_work<int>(8,5);
    mult_work<long>(8,5);
//
//
//
//    mult_work<double>(5,8);
//    mult_work<float>(5,8);
//    mult_work<int>(5,8);
//    mult_work<long>(5,8);


}
