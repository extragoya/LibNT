#include <string>

#define BOOST_TEST_MODULE LatticeConstrTests



#include "DenseLattice.h"
#include "SparseLattice.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

constexpr int dim=5;

template<typename T>
void constructor_work(){
    typedef T Lat1;


    Lat1 test;
    BOOST_CHECK_EQUAL(test.height(),0);
    BOOST_CHECK_EQUAL(test.width(),0);
    BOOST_CHECK_EQUAL(test.depth(),0);
    BOOST_CHECK_EQUAL(test.data_end()-test.data_begin(),0);

    Lat1 test2(dim,dim,dim);
    BOOST_CHECK_EQUAL(test2.height(),dim);
    BOOST_CHECK_EQUAL(test2.width(),dim);
    BOOST_CHECK_EQUAL(test2.depth(),dim);





}

template<typename T>
void dense_constructor_work(){
    typedef T Lat1;
    constructor_work<T>();
    Lat1 test2(dim,dim,dim);

    BOOST_CHECK_EQUAL(test2.data_end()-test2.data_begin(),dim*dim*dim);
    for(auto it=test2.data_begin();it<test2.data_end();++it)
        BOOST_CHECK_EQUAL(*it,0);

    test2.randu(-5,5);
    Lat1 test3(test2.data_begin(),test2.height(),test2.width(),test2.depth());
    BOOST_CHECK_EQUAL(test3.height(),dim);
    BOOST_CHECK_EQUAL(test3.width(),dim);
    BOOST_CHECK_EQUAL(test3.depth(),dim);
    for(auto it=test2.data_begin(),it2=test3.data_begin();it<test2.data_end();++it,++it2)
        BOOST_CHECK_EQUAL(*it,*it2);

}

template<typename T>
void sparse_constructor_work(){
    typedef T Lat1;

    constructor_work<T>();

    Lat1 test;
    BOOST_CHECK_EQUAL(test.index_end()-test.index_begin(),0);

    Lat1 test2(dim,dim,dim);
    BOOST_CHECK_EQUAL(test2.data_end()-test2.data_begin(),0);
    BOOST_CHECK_EQUAL(test2.index_end()-test2.index_begin(),0);
    test2.eye();
    Lat1 test3(test2.data(),test2.indices(),test2.height(),test2.width(),test2.depth());




}

BOOST_AUTO_TEST_CASE( DenseLatticeConstrTests )
{

    dense_constructor_work<LibMIA::DenseLattice<float> >();
    dense_constructor_work<LibMIA::DenseLattice<double> >();
    dense_constructor_work<LibMIA::DenseLattice<int32_t> >();
    dense_constructor_work<LibMIA::DenseLattice<int64_t> >();
}


BOOST_AUTO_TEST_CASE( SparseLatticeConstrTests )
{

    sparse_constructor_work<LibMIA::SparseLattice<float> >();
    sparse_constructor_work<LibMIA::SparseLattice<double> >();
    sparse_constructor_work<LibMIA::SparseLattice<int32_t> >();
    sparse_constructor_work<LibMIA::SparseLattice<int64_t> >();
}
