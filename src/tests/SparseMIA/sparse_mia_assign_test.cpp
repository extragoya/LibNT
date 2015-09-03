#include <iostream>

#define BOOST_TEST_MODULE SparseAssignTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "SparseMIA.h"

#include "Index.h"

//typedef LibMIA::DenseMIA<double,3> dmia;

template<class MIAType>
void assign_work(){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;

    MIAType a(5,7,4);
    MIAType b;
    MIAType c;


    c=a;
    BOOST_CHECK_MESSAGE(a==c && c==a,"Empty SparseMIA assignment");

    a.resize(a.dimensionality()/2);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    c=a;
    BOOST_CHECK_MESSAGE(a==c && c==a,"Straight SparseMIA assignment");

    b.zeros();
    c.zeros();
    b(i,j,k)=a(i,j,k);
    c(i,j,k)=b(i,j,k);
    BOOST_CHECK_MESSAGE(a==c,"Non-shuffled MIA Expression assignment");

    b.zeros();
    c.zeros();
    b(i,j,k)=a(i,k,j);
    c(i,j,k)=b(i,k,j);
    BOOST_CHECK_MESSAGE(a==c,"Shuffled MIA Expression assignment");

    b.zeros();
    c.zeros();
    b(i,j,k)=a(i,!k,j);
    c(!i,j,k)=b(i,k,j);
    BOOST_CHECK_MESSAGE(a==c,"Shuffled MIA Expression assignment with mismatched elem-wise indices");

}

BOOST_AUTO_TEST_CASE( SparseAssignTests )
{

    assign_work<LibMIA::SparseMIA<double,3>>();
    assign_work<LibMIA::SparseMIA<float,3>>();
    assign_work<LibMIA::SparseMIA<int,3>>();
    assign_work<LibMIA::SparseMIA<long,3>>();



}
