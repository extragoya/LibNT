#include <iostream>

#define BOOST_TEST_MODULE MIAAssignTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"

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

    a.randu(-5,5);
    const MIAType temp_a(a); //check const correctness
    c=temp_a;
    BOOST_CHECK_MESSAGE(a==c,"Straight MIA assignment");

    b.zeros();
    c.zeros();
    b(i,j,k)=temp_a(i,j,k);
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

BOOST_AUTO_TEST_CASE( MIAAssignTests )
{

    assign_work<LibMIA::DenseMIA<double,3>>();
    assign_work<LibMIA::DenseMIA<float,3>>();
    assign_work<LibMIA::DenseMIA<int,3>>();
    assign_work<LibMIA::DenseMIA<long,3>>();



}
