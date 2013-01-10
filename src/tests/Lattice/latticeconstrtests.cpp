#include <string>

#define BOOST_TEST_MODULE DenseLatticeConstrTests
#include <boost/test/included/unit_test.hpp>


#include "DenseLattice.h"



template<typename T1>
void constructor_work(){
    typedef LibMIA::DenseLattice<T1> Lat1;


    Lat1 test;
    BOOST_CHECK_EQUAL(test.height(),0);
    BOOST_CHECK_EQUAL(test.width(),0);
    BOOST_CHECK_EQUAL(test.depth(),0);
    BOOST_CHECK_EQUAL(test.data_end()-test.data_begin(),0);



}

BOOST_AUTO_TEST_CASE( DenseLatticeConstrTests )
{
    constructor_work<float>();
    constructor_work<double>();
    constructor_work<int32_t>();
    constructor_work<int64_t>();
}
