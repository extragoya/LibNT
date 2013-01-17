#include <iostream>

#define BOOST_TEST_MODULE MIAAssignTests
#include <boost/test/included/unit_test.hpp>

#include "DenseMIA.h"
#include "Index.h"



typedef LibMIA::DenseMIA<double,3> dmia;

BOOST_AUTO_TEST_CASE( MIAAssignTests )
{

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;

    dmia a(5,7,4);
    dmia b;
    dmia c;
    a.randu(-5,5);

    c=a;
    BOOST_CHECK_MESSAGE(a==c,"Straight MIA assignment");


    b(i,j,k)=a(i,j,k);
    c(i,j,k)=b(i,j,k);
    BOOST_CHECK_MESSAGE(a==c,"Non-shuffled MIA Expression assignment");

    b(i,j,k)=a(i,k,j);
    c(i,j,k)=b(i,k,j);
    BOOST_CHECK_MESSAGE(a==c,"Shuffled MIA Expression assignment");


}
