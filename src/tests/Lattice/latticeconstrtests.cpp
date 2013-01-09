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

    /*//= Lat1(5,5,10);
    test1.load("data/test1.bin");
    Lat2 test2= Lat2(5,5,10);

    test2.load("data/test2.bin");
    Lat3 test3=test1*test2;
    test3.save(_result);
    Lat3 test3check;
    test3check.load(_result);
    BOOST_CHECK( test3 == test3check );*/

}

BOOST_AUTO_TEST_CASE( DenseLatticeConstrTests )
{


    constructor_work<float>();
    constructor_work<double>();
    constructor_work<int32_t>();
    constructor_work<int64_t>();





}
