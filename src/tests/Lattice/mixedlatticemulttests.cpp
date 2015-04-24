#define BOOST_TEST_MODULE MixedLatticeMultTests



#include "DenseLattice.h"
#include "SparseLattice.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif


template<typename data_type>
void multwork(size_t m1, size_t n1, size_t n2, size_t p){
    typedef LibMIA::DenseLattice<data_type> denseType;
    typedef LibMIA::SparseLattice<data_type> sparseType;
    denseType DenseLat1(m1,n1,p);
    denseType DenseLat2(n1,n2,p);

    denseType DenseLat3_dense;
    sparseType DenseLat3_sparse;

    sparseType SparseLat1(m1,n1,p);
    sparseType SparseLat2(n1,n2,p);




    DenseLat1.randu(0,10);
    DenseLat2.randu(0,10);
    for(auto it=DenseLat1.data_begin();it<DenseLat1.data_end();++it)
        if(*it<7)
            *it=0;
    for(auto it=DenseLat2.data_begin();it<DenseLat2.data_end();++it)
        if(*it<7)
            *it=0;


    SparseLat1=DenseLat1;
    SparseLat2=DenseLat2;


    DenseLat3_dense=DenseLat1*DenseLat2;

    DenseLat3_sparse=SparseLat1*DenseLat2;


//    DenseLat3_dense.print();
//    DenseLat3_sparse.print();
    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Full Dimension Mult Test 1 for ")+typeid(data_type).name());


    DenseLat3_sparse=DenseLat1*SparseLat2;

    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Full Dimension Mult Test 2 for ")+typeid(data_type).name());

    DenseLat1=denseType(1,n1,p);
    DenseLat2=denseType(n1,1,p);
    DenseLat1.randu(0,10);
    DenseLat2.randu(0,10);
    for(auto it=DenseLat1.data_begin();it<DenseLat1.data_end();++it)
        if(*it<7)
            *it=0;
    for(auto it=DenseLat2.data_begin();it<DenseLat2.data_end();++it)
        if(*it<7)
            *it=0;

    SparseLat1=DenseLat1;
    SparseLat2=DenseLat2;
    //std::cout << "About to dense mult" << std::endl;
    DenseLat3_dense=DenseLat1*DenseLat2;
    //std::cout << "About to mixed mult" << std::endl;
    DenseLat3_sparse=SparseLat1*DenseLat2;

    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Repeated Inner Product Mult Test 1 for ")+typeid(data_type).name());

    DenseLat3_sparse=DenseLat1*SparseLat2;

    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Repeated Inner Product Mult Test 2 for ")+typeid(data_type).name());

    DenseLat1=denseType(m1,1,p);
    DenseLat2=denseType(1,n2,p);
    DenseLat1.randu(0,10);
    DenseLat2.randu(0,10);
    for(auto it=DenseLat1.data_begin();it<DenseLat1.data_end();++it)
        if(*it<7)
            *it=0;
    for(auto it=DenseLat2.data_begin();it<DenseLat2.data_end();++it)
        if(*it<7)
            *it=0;

    SparseLat1=DenseLat1;
    SparseLat2=DenseLat2;

    DenseLat3_dense=DenseLat1*DenseLat2;
    DenseLat3_sparse=SparseLat1*DenseLat2;

    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Repeated Outer Product Mult Test 1 for ")+typeid(data_type).name());

    DenseLat3_sparse=DenseLat1*SparseLat2;

    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Repeated Outer Product Mult Test 2 for ")+typeid(data_type).name());

    DenseLat1=denseType(m1,n1,1);
    DenseLat2=denseType(n1,n2,1);
    DenseLat1.randu(0,10);
    DenseLat2.randu(0,10);
    for(auto it=DenseLat1.data_begin();it<DenseLat1.data_end();++it)
        if(*it<7)
            *it=0;
    for(auto it=DenseLat2.data_begin();it<DenseLat2.data_end();++it)
        if(*it<7)
            *it=0;

    SparseLat1=DenseLat1;
    SparseLat2=DenseLat2;

    DenseLat3_dense=DenseLat1*DenseLat2;
    DenseLat3_sparse=SparseLat1*DenseLat2;
    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("No Depth Mult Test 1 for ")+typeid(data_type).name());
    DenseLat3_sparse=DenseLat1*SparseLat2;
    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("No Depth Mult Test 2 for ")+typeid(data_type).name());



    DenseLat1=denseType(m1,1,p);
    DenseLat2=denseType(1,1,p);
    DenseLat1.randu(0,10);
    DenseLat2.randu(0,10);

    for(auto it=DenseLat1.data_begin();it<DenseLat1.data_end();++it)
        if(*it<7)
            *it=0;
    for(auto it=DenseLat2.data_begin();it<DenseLat2.data_end();++it)
        if(*it<7)
            *it=0;

    SparseLat1=DenseLat1;
    SparseLat2=DenseLat2;

    DenseLat3_dense=DenseLat1*DenseLat2;
    DenseLat3_sparse=SparseLat1*DenseLat2;

    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Outer product, one operand, test 1 for ")+typeid(data_type).name());
    DenseLat3_sparse=DenseLat1*SparseLat2;
    BOOST_CHECK_MESSAGE(DenseLat3_dense.fuzzy_equals(DenseLat3_sparse,test_precision<data_type>()),std::string("Outer product, one operand, test 2 for ")+typeid(data_type).name());




}

BOOST_AUTO_TEST_CASE( SparseLatticeMultTests )
{


    //multwork<double>(5,5,5,5);
    multwork<double>(20,20,20,20);
    multwork<float>(20,20,20,20);
    multwork<int>(20,20,20,20);
    multwork<long>(20,20,20,20);



}
