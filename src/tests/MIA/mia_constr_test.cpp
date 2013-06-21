#include <string>

#define BOOST_TEST_MODULE MIAConstrTests




#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"
#include "SparseMIA.h"
constexpr int dim=5;



template<typename MIAType>
void constructor_work(size_t dim1, size_t dim2, size_t dim3){



    MIAType test;
    BOOST_CHECK_MESSAGE(test.dims()[0]==0,std::string("Zero-size dimension test 1 for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test.dims()[1]==0,std::string("Zero-size dimension test 2 for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test.dims()[2]==0,std::string("Zero-size dimension test 3 for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test.dimensionality()==0,std::string("Zero dimensionality test for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test.data_end()-test.data_begin()==0,std::string("Zero-size scalar size for ")+typeid(MIAType).name());

    MIAType test2(dim1,dim2,dim3);
    BOOST_CHECK_MESSAGE(test2.dims()[0]==dim1,std::string("Size dimension test 1 for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test2.dims()[1]==dim2,std::string("Size dimension test 1 for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test2.dims()[2]==dim3,std::string("Size dimension test 1 for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test2.dimensionality()==dim1*dim2*dim3,std::string("Full dimensionality test for ")+typeid(MIAType).name());





}

template<typename MIAType>
void dense_constructor_work(size_t dim1, size_t dim2, size_t dim3){

    constructor_work<MIAType>(dim1,dim2,dim3);
    MIAType test2(dim1,dim2,dim3);

    BOOST_CHECK_EQUAL(test2.data_end()-test2.data_begin(),(int)(dim1*dim2*dim3));
    bool passed=true;
    for(auto it=test2.data_begin();it<test2.data_end();++it)
        passed=(passed&&*it==0);

    BOOST_CHECK_MESSAGE(passed,"Constructor initialized to zero");

    test2.randu(-5,5);
    MIAType test3(test2.dims(),test2.raw_data_ptr(),true);
    test2.release_raw_data();
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

template<typename MIAType, typename DenseMIAType>
void sparse_constructor_work(size_t dim1, size_t dim2, size_t dim3){


    constructor_work<MIAType>(dim1,dim2,dim3);
    //empty constructor - test index size (data size is tested in constructor work)
    MIAType test1;
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==0,std::string("Zero-size index size with no dimensionality for ")+typeid(MIAType).name());

    //test size of data and index containers when given dimensionality - should be zero as all values are zero
    MIAType test2(dim1,dim2,dim3);
    BOOST_CHECK_MESSAGE(test2.data_end()-test2.data_begin()==0,std::string("Zero-size data size with dimensionality for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test2.index_end()-test2.index_begin()==0,std::string("Zero-size index size with dimensionality for ")+typeid(MIAType).name());

    //create a densemia for copy constructor tests
    DenseMIAType denseTest(dim1,dim2,dim3);
    //create one with no zeros
    denseTest.randu(1,2);
    MIAType test3(denseTest);
    //check that data and index containers are completely full
    BOOST_CHECK_MESSAGE(test3.index_end()-test3.index_begin()==(int)(dim1*dim2*dim3),std::string("Dimensionality-size index test for DenseMIA with no zeros copy constructor ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test3.data_end()-test3.data_begin()==(int)(dim1*dim2*dim3),std::string("Dimensionality-size data test for DenseMIA with no zeros copy constructor ")+typeid(MIAType).name());
    //now check that the contents match that of the denseMIA
    bool passed=true;
    for(auto it=denseTest.data_begin(),it2=test3.data_begin();it<denseTest.data_end();++it,++it2)
        passed=(passed&&*it==*it2);
    BOOST_CHECK_MESSAGE(passed,std::string("Data content test for DenseMIA with no zeros copy constructor for ")+typeid(MIAType).name());
    //now test that the index is a sequence of integers going from 0 to dimensionality-1
    passed=true;
    for(auto it2=test3.index_begin();it2<test3.index_end();++it2)
        passed=(passed&&*it2==(it2-test3.index_begin()));
    BOOST_CHECK_MESSAGE(passed,std::string("Index content test for DenseMIA with no zeros copy constructor for ")+typeid(MIAType).name());


    denseTest.randu(-1,1);
    for(auto it=denseTest.data_begin();it<denseTest.data_end();++it)
        if(*it<0)
            *it=0;
    MIAType test4(denseTest);
    BOOST_CHECK_MESSAGE(test4==denseTest,std::string("Sparse copy constructor for DenseMIA for ")+typeid(MIAType).name());
    denseTest.zeros();
    MIAType test5(denseTest);
    passed=test5.index_end()-test5.index_begin()==0;
    passed=passed&&test5.data_end()-test5.data_begin()==0;
    passed=passed&&test5.size()==0;
    BOOST_CHECK_MESSAGE(passed,std::string("Sparse copy constructor for all-zero DenseMIA for ")+typeid(MIAType).name());

    test5.resize(test5.dimensionality()/2);
    test5.randu(-50,50);
    test5.rand_indices();
    test5.collect_duplicates();

    MIAType test6(test5);
    passed=test6.linIdxSequence()==test5.linIdxSequence();
    passed=passed&&test5.index_end()-test5.index_begin()==test6.index_end()-test6.index_begin();
    passed=passed&&test5.data_end()-test5.data_begin()==test6.data_end()-test6.data_begin();
    passed=passed&&test5.size()==test6.size();
    auto it5=test5.data_begin();
    for(auto it=test6.data_begin();it<test6.data_end();++it){
        passed=passed&&*it==*it5;
        it5++;
    }
    auto it5_2=test5.index_begin();
    for(auto it=test6.index_begin();it<test6.index_end();++it){
        passed=passed&&*it==*it5_2;
        it5_2++;
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Sparse copy constructor for SparseMIA for ")+typeid(MIAType).name());



}

BOOST_AUTO_TEST_CASE( DenseMIAConstrTests )
{

    size_t dim1=3,dim2=4,dim3=5;


    dense_constructor_work<LibMIA::DenseMIA<float,3> >(dim1,dim2,dim3);
    dense_constructor_work<LibMIA::DenseMIA<double,3> >(dim1,dim2,dim3);
    dense_constructor_work<LibMIA::DenseMIA<int32_t,3> >(dim1,dim2,dim3);
    dense_constructor_work<LibMIA::DenseMIA<int64_t,3> >(dim1,dim2,dim3);

}


BOOST_AUTO_TEST_CASE( SparseMIAConstrTests )
{

    size_t dim1=3,dim2=4,dim3=5;
    //testing to make sure a copy constructor doesn't get called erroneously
    LibMIA::SparseMIA<float,1> test(dim1);
    sparse_constructor_work<LibMIA::SparseMIA<float,3>, LibMIA::DenseMIA<float,3> >(dim1,dim2,dim3);
    sparse_constructor_work<LibMIA::SparseMIA<double,3>, LibMIA::DenseMIA<double,3> >(dim1,dim2,dim3);
    sparse_constructor_work<LibMIA::SparseMIA<int32_t,3>,LibMIA::DenseMIA<int32_t,3> >(dim1,dim2,dim3);
    sparse_constructor_work<LibMIA::SparseMIA<int64_t,3>,LibMIA::DenseMIA<int64_t,3> >(dim1,dim2,dim3);
}
