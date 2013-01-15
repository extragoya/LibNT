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
    dmia b(7,4,5);
    dmia c(5,7,4);

    a.randu(-5,5);
    b(i,j,k)=a(k,i,j);
    c(i,j,k)=b(j,k,i);
    BOOST_CHECK(a==b);


}
