#include <string>

#define BOOST_TEST_MODULE MIAConstrTests



#include "DenseMIA.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

constexpr int dim=5;

template<typename MIAType>
void constructor_work(size_t dim1, size_t dim2, size_t dim3){



    MIAType test;
    BOOST_CHECK_EQUAL(test.dims()[0],0);
    BOOST_CHECK_EQUAL(test.dims()[1],0);
    BOOST_CHECK_EQUAL(test.dims()[2],0);
    BOOST_CHECK_EQUAL(test.data_end()-test.data_begin(),0);

    MIAType test2(dim1,dim2,dim3);
    BOOST_CHECK_EQUAL(test2.dims()[0],dim1);
    BOOST_CHECK_EQUAL(test2.dims()[1],dim2);
    BOOST_CHECK_EQUAL(test2.dims()[2],dim3);





}

template<typename MIAType>
void dense_constructor_work(size_t dim1, size_t dim2, size_t dim3){

    constructor_work<MIAType>(dim1,dim2,dim3);
    MIAType test2(dim1,dim2,dim3);

    BOOST_CHECK_EQUAL(test2.data_end()-test2.data_begin(),dim1*dim2*dim3);
    bool passed=true;
    for(auto it=test2.data_begin();it<test2.data_end();++it)
        passed=(passed&&*it==0);

    BOOST_CHECK_MESSAGE(passed,"Constructor initialized to zero");

    test2.randu(-5,5);
    MIAType test3(test2.dims(),test2.raw_data_ptr(),true);
    test2.release_ownership();
    BOOST_CHECK_EQUAL(test3.dims()[0],dim1);
    BOOST_CHECK_EQUAL(test3.dims()[1],dim2);
    BOOST_CHECK_EQUAL(test3.dims()[2],dim3);
    passed=true;
    for(auto it=test2.data_begin(),it2=test3.data_begin();it<test2.data_end();++it,++it2)
        passed=(passed&&it==it2);

    BOOST_CHECK_MESSAGE(passed,"Constructor passing pointer without ownership");

    MIAType test4(test2.dims(),test2.raw_data_ptr(),false);
    BOOST_CHECK_EQUAL(test4.dims()[0],dim1);
    BOOST_CHECK_EQUAL(test4.dims()[1],dim2);
    BOOST_CHECK_EQUAL(test4.dims()[2],dim3);
    passed=true;
    for(auto it=test2.data_begin(),it2=test4.data_begin();it<test2.data_end();++it,++it2)
        passed=(passed&&*it==*it2);

    BOOST_CHECK_MESSAGE(passed,"Constructor passing pointer without ownership");

}

//template<typename T>
//void sparse_constructor_work(){
//    typedef T Lat1;
//
//    constructor_work<T>();
//
//    Lat1 test;
//    BOOST_CHECK_EQUAL(test.index_end()-test.index_begin(),0);
//
//    Lat1 test2(dim,dim,dim);
//    BOOST_CHECK_EQUAL(test2.data_end()-test2.data_begin(),0);
//    BOOST_CHECK_EQUAL(test2.index_end()-test2.index_begin(),0);
//    test2.eye();
//    Lat1 test3(test2.data(),test2.indices(),test2.height(),test2.width(),test2.depth());
//
//
//
//
//}

BOOST_AUTO_TEST_CASE( DenseMIAConstrTests )
{

    size_t dim1=3,dim2=4,dim3=5;


    dense_constructor_work<LibMIA::DenseMIA<float,3> >(dim1,dim2,dim3);
    dense_constructor_work<LibMIA::DenseMIA<double,3> >(dim1,dim2,dim3);
    dense_constructor_work<LibMIA::DenseMIA<int32_t,3> >(dim1,dim2,dim3);
    dense_constructor_work<LibMIA::DenseMIA<int64_t,3> >(dim1,dim2,dim3);

}


//BOOST_AUTO_TEST_CASE( SparseLatticeConstrTests )
//{
//
//    sparse_constructor_work<LibMIA::SparseLattice<float> >();
//    sparse_constructor_work<LibMIA::SparseLattice<double> >();
//    sparse_constructor_work<LibMIA::SparseLattice<int32_t> >();
//    sparse_constructor_work<LibMIA::SparseLattice<int64_t> >();
//}
