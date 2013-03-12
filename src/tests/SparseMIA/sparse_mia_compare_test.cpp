#include <string>

#define BOOST_TEST_MODULE SparseCompareTests




#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"
#include "SparseMIA.h"






template<typename MIAType, typename DenseMIAType>
void do_work(size_t dim1, size_t dim2, size_t dim3){


    MIAType test1(dim1,dim2,dim3), test2(dim1,dim2,dim3);
    test1.resize(2);
    test2.resize(2);
    std::fill(test1.data_begin(),test1.data_end(),1);
    std::fill(test2.data_begin(),test2.data_end(),1);
    *(test1.index_begin())=1;*(test2.index_begin())=1;
    *(test1.index_begin()+1)=2;*(test2.index_begin()+1)=2;
    BOOST_CHECK_MESSAGE(test1==test2,std::string("Comparison test 1 for ")+typeid(MIAType).name());
    *(test1.index_begin()+1)=2;*(test2.index_begin()+1)=3;
    BOOST_CHECK_MESSAGE(test1!=test2,std::string("Comparison test 2 for ")+typeid(MIAType).name());
    *(test1.index_begin()+1)=2;*(test2.index_begin()+1)=2;
    *(test1.data_begin())=2;
    BOOST_CHECK_MESSAGE(test1!=test2,std::string("Comparison test 3 for ")+typeid(MIAType).name());

    test1.resize(test1.dimensionality()/2);
    test2.resize(test1.dimensionality()/2);
    test1.randu(-50,50);
    test1.rand_indices();
    test1.collect_duplicates();

    std::copy(test1.storage_begin(),test1.storage_end(),test2.storage_begin());
    BOOST_CHECK_MESSAGE(test1==test2,std::string("Comparison test 4 for ")+typeid(MIAType).name());




}




BOOST_AUTO_TEST_CASE( SparseCompareTests )
{

    size_t dim1=3,dim2=4,dim3=5;
    //testing to make sure a copy constructor doesn't get called erroneously
    do_work<LibMIA::SparseMIA<float,3>, LibMIA::DenseMIA<float,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<double,3>, LibMIA::DenseMIA<double,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<int32_t,3>,LibMIA::DenseMIA<int32_t,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<int64_t,3>,LibMIA::DenseMIA<int64_t,3> >(dim1,dim2,dim3);
}
