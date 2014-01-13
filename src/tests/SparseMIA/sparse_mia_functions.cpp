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






template<typename data_type>
void do_work(size_t dim1, size_t dim2, size_t dim3){


    typedef LibMIA::SparseMIA<data_type,3> MIAType;
    MIAType test1(dim1,dim2,dim3);
    test1.resize(test1.dimensionality()/2);
    //testing resize functionality
    BOOST_CHECK_MESSAGE(test1.data_end()-test1.data_begin()==test1.dimensionality()/2,std::string("Resize test of data for ")+typeid(data_type).name());
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==test1.dimensionality()/2,std::string("Resize test of indices for ")+typeid(data_type).name());

    //test sorting
    test1.rand_indices();
    test1.sort();
    bool passed=true;
    for(auto i=test1.index_begin()+1;i<test1.index_end();++i)
        passed=(passed&& *(i-1)<=*i);
    BOOST_CHECK_MESSAGE(passed,std::string("Sorting test for ")+typeid(data_type).name());

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
    BOOST_CHECK_MESSAGE(test1.data_end()-test1.data_begin()==2,std::string("Collect duplicates test 1 of data size for ")+typeid(data_type).name());
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==2,std::string("Collect duplicates test 1 of indices size for ")+typeid(data_type).name());
    passed=(*(test1.data_begin())==1 && *(test1.data_begin()+1)==2);
    BOOST_CHECK_MESSAGE(passed,std::string("Collect duplicates test 1 of data content for ")+typeid(data_type).name());

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


    test1.collect_duplicates(std::plus<data_type>());
    BOOST_CHECK_MESSAGE(test1.data_end()-test1.data_begin()==2,std::string("Collect duplicates test 2 of data size for ")+typeid(data_type).name());
    BOOST_CHECK_MESSAGE(test1.index_end()-test1.index_begin()==2,std::string("Collect duplicates test 2 of indices size for ")+typeid(data_type).name());
    passed=(*(test1.data_begin())==(data_type)counter1 && *(test1.data_begin()+1)==(data_type)counter2);
    BOOST_CHECK_MESSAGE(passed,std::string("Collect duplicates test 2 of data content for ")+typeid(data_type).name());

    MIAType test2; //zero dimensionality

    test2.rand_indices(); //make sure no error thrown when dimensionality is zero


    LibMIA::SparseMIA<data_type,4> test_reorder(dim1,dim2,dim3,dim2);
    LibMIA::SparseMIA<data_type,4> test_reorder2;
    auto linIdxSequence=test_reorder.linIdxSequence();

    linIdxSequence={{2,0,3,1}};

    test_reorder.resize(test_reorder.dimensionality()/2);
    test_reorder.randu(-5,5);
    test_reorder.rand_indices();
    test_reorder.collect_duplicates();
    test_reorder2=test_reorder;
    test_reorder.sort(linIdxSequence);
    test_reorder2.change_linIdx_sequence(linIdxSequence);
    test_reorder2.sort();
    passed=true;
    auto data_it=test_reorder.data_begin();
    auto data_it2=test_reorder2.data_begin();
    auto it=test_reorder.index_begin();
    auto it2=test_reorder2.index_begin();
    for(;data_it<test_reorder.data_end();++data_it,++it,++data_it2,++it2){ //do a manual check, b/c if we use == operator, a sort will occur, nullifying the test
        if(*data_it!=*data_it2){
            passed=false;
            break;
        }
        if(*it!=*it2){
            passed=false;
            break;
        }
    }

    BOOST_CHECK_MESSAGE(passed,std::string("Sparse reording sort 1 for ")+typeid(data_type).name());

    //try another linIdxSequence
    linIdxSequence={{1,2,3,0}};
    test_reorder2=test_reorder;
    test_reorder.sort(linIdxSequence);
    test_reorder2.change_linIdx_sequence(linIdxSequence);
    test_reorder2.sort();
    passed=true;
    data_it=test_reorder.data_begin();
    data_it2=test_reorder2.data_begin();
    it=test_reorder.index_begin();
    it2=test_reorder2.index_begin();
    for(;data_it<test_reorder.data_end();++data_it,++it,++data_it2,++it2){ //do a manual check, b/c if we use == operator, a sort will occur, nullifying the test
        if(*data_it!=*data_it2){
            passed=false;
            break;
        }
        if(*it!=*it2){
            passed=false;
            break;
        }
    }

    BOOST_CHECK_MESSAGE(passed,std::string("Sparse reording sort 2 for ")+typeid(data_type).name());

    //try another linIdxSequence
    linIdxSequence={{0,1,2,3}};
    test_reorder2=test_reorder;
    test_reorder.sort(linIdxSequence);
    test_reorder2.change_linIdx_sequence(linIdxSequence);
    test_reorder2.sort();
    passed=true;
    data_it=test_reorder.data_begin();
    data_it2=test_reorder2.data_begin();
    it=test_reorder.index_begin();
    it2=test_reorder2.index_begin();
    for(;data_it<test_reorder.data_end();++data_it,++it,++data_it2,++it2){ //do a manual check, b/c if we use == operator, a sort will occur, nullifying the test
        if(*data_it!=*data_it2){
            passed=false;
            break;
        }
        if(*it!=*it2){
            passed=false;
            break;
        }
    }

    BOOST_CHECK_MESSAGE(passed,std::string("Sparse reording sort 3 for ")+typeid(data_type).name());

    //try another linIdxSequence
    linIdxSequence={{2,3,0,1}};
    test_reorder2=test_reorder;
    test_reorder.sort(linIdxSequence);
    test_reorder2.change_linIdx_sequence(linIdxSequence);
    test_reorder2.sort();
    passed=true;
    data_it=test_reorder.data_begin();
    data_it2=test_reorder2.data_begin();
    it=test_reorder.index_begin();
    it2=test_reorder2.index_begin();
    for(;data_it<test_reorder.data_end();++data_it,++it,++data_it2,++it2){ //do a manual check, b/c if we use == operator, a sort will occur, nullifying the test
        if(*data_it!=*data_it2){
            passed=false;
            break;
        }
        if(*it!=*it2){
            passed=false;
            break;
        }
    }

    BOOST_CHECK_MESSAGE(passed,std::string("Sparse reording sort 3 for ")+typeid(data_type).name());


}




BOOST_AUTO_TEST_CASE( SparseCompareTests )
{

    //size_t dim1=3,dim2=4,dim3=5;
    size_t dim1=3,dim2=3,dim3=3;


    do_work<float>(dim1,dim2,dim3);
//    do_work<double>(dim1,dim2,dim3);
//    do_work<int32_t>(dim1,dim2,dim3);
//    do_work<int64_t>(dim1,dim2,dim3);
}
