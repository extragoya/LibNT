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






template<typename MIAType>
void do_work(size_t dim1, size_t dim2, size_t dim3){


    MIAType test1(dim1,dim2,dim3);
    test1.resize(test1.dimensionality()/2);
    //testing resize functionality
    BOOST_CHECK_MESSAGE(test1.data_end()-test1.data_begin()==test1.dimensionality()/2,std::string("Resize test of data for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==test1.dimensionality()/2,std::string("Resize test of indices for ")+typeid(MIAType).name());

    //test sorting
    test1.rand_indices();
    test1.sort();
    bool passed=true;
    for(auto i=test1.index_begin()+1;i<test1.index_end();++i)
        passed=(passed&& *(i-1)<=*i);
    BOOST_CHECK_MESSAGE(passed,std::string("Sorting test for ")+typeid(MIAType).name());

    //testing removal of duplicates
    test1.randu(0,2);
    *(test1.data_begin())=2;
    *(test1.data_begin()+1)=1;
    *(test1.index_begin())=2;
    *(test1.index_begin()+1)=1;
    test1.setSorted(false);
    auto i=test1.index_begin()+2;
    for(auto data_it=test1.data_begin()+2;i<test1.index_end();++i,++data_it)
        if (*data_it>=1){
            *data_it=2;
            *i=2;
        }
        else{
            *data_it=1;
            *i=1;
        }

    test1.collect_duplicates();
    BOOST_CHECK_MESSAGE(test1.data_end()-test1.data_begin()==2,std::string("Collect duplicates test 1 of data size for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==2,std::string("Collect duplicates test 1 of indices size for ")+typeid(MIAType).name());
    passed=(*(test1.data_begin())==1 && *(test1.data_begin()+1)==2);
    BOOST_CHECK_MESSAGE(passed,std::string("Collect duplicates test 1 of data content for ")+typeid(MIAType).name());

    //testing collection of duplicates where we add duplicated entries
    test1.resize(test1.dimensionality()/2);
    test1.randu(0,2);
    *(test1.data_begin())=1;
    *(test1.data_begin()+1)=1;
    *(test1.index_begin())=2;
    *(test1.index_begin()+1)=1;
    test1.setSorted(false);
    size_t counter1=1,counter2=1;
    i=test1.index_begin()+2;
    for(auto data_it=test1.data_begin()+2;i<test1.index_end();++i,++data_it)
        if (*data_it>=1){
            *data_it=1;
            *i=1;
            ++counter1;
        }
        else{
            *data_it=1;
            *i=2;
            ++counter2;
        }

    typedef typename MIAType::data_type data_type;
    test1.collect_duplicates(std::plus<data_type>());
    BOOST_CHECK_MESSAGE(test1.data_end()-test1.data_begin()==2,std::string("Collect duplicates test 2 of data size for ")+typeid(MIAType).name());
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==2,std::string("Collect duplicates test 2 of indices size for ")+typeid(MIAType).name());
    passed=(*(test1.data_begin())==(data_type)counter1 && *(test1.data_begin()+1)==(data_type)counter2);
    BOOST_CHECK_MESSAGE(passed,std::string("Collect duplicates test 2 of data content for ")+typeid(MIAType).name());

    MIAType test2; //zero dimensionality

    test2.rand_indices(); //make sure no error thrown when dimensionality is zero

}




BOOST_AUTO_TEST_CASE( SparseCompareTests )
{

    size_t dim1=3,dim2=4,dim3=5;


    do_work<LibMIA::SparseMIA<float,3>>(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<double,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<int32_t,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<int64_t,3> >(dim1,dim2,dim3);
}
