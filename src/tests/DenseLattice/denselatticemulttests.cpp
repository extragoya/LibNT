#include <string>

#define BOOST_TEST_MODULE DenseLatticeMultTests



#include "DenseLattice.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif


template<typename T1,typename T2,typename T3>
void multwork(const std::string & _result){
    typedef LibMIA::DenseLattice<T1> Lat1;
    typedef LibMIA::DenseLattice<T2> Lat2;
    typedef LibMIA::DenseLattice<T3> Lat3;

    Lat1 test1= Lat1(5,5,10);
    test1.load("data/test1.bin");
    Lat2 test2= Lat2(5,5,10);

    test2.load("data/test2.bin");
    Lat3 test3=test1*test2;
    test3.save(_result);
    Lat3 test3check;
    test3check.load(_result);
    BOOST_CHECK( test3 == test3check );

}

BOOST_AUTO_TEST_CASE( DenseLatticeMultTests )
{


    multwork<double,double,double>("data/multtest_ddd.bin");
    multwork<float,float,double>("data/multtest_dff.bin");
    multwork<float,float,float>("data/multtest_fff.bin");
    multwork<double,double,float>("data/multtest_ddf.bin");
    multwork<int32_t,int32_t,int32_t>("data/multtest_iii.bin");
    multwork<int64_t,int64_t,int64_t>("data/multtest_lll.bin");




}
