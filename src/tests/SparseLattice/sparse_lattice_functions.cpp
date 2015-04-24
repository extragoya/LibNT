#define BOOST_TEST_MODULE SparseLatticeFunctionTests




#include "SparseLattice.h"
#include "DenseLattice.h"

#include "MIAConfig.h"



#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
//#include <boost/timer/timer.hpp>

template<typename data_type>
void do_work(size_t m, size_t n, size_t p){

    typedef LibMIA::SparseLattice<data_type> sparseType;
    typedef LibMIA::DenseLattice<data_type> denseType;
    denseType denseLat(m,n,p);
    sparseType sparseLat1,sparseLat2;

    denseLat.randu(0,10);
    for(auto it=denseLat.data_begin();it<denseLat.data_end();++it)
        if(*it<7)
            *it=0;
    sparseLat1=denseLat;
    sparseLat2=denseLat;

    BOOST_CHECK_MESSAGE(sparseLat1==sparseLat2,std::string("Assignment and comparison test for ")+typeid(data_type).name());


    //boost::timer::cpu_timer ltensor;
    sparseLat1.sort(LibMIA::RowMajor);
    //std::cout << "special sort time " << boost::timer::format(ltensor.elapsed()) << std::endl;




    sparseLat2.set_sorted(false);
    //ltensor=boost::timer::cpu_timer ();

    sparseLat2.sort(LibMIA::RowMajor);
    //std::cout << "normal sort time " << boost::timer::format(ltensor.elapsed()) << std::endl;
    //sparseLat1.print();
    //sparseLat2.print();

    //sparseLat1.print();
    //sparseLat2.print();

    BOOST_CHECK_MESSAGE(sparseLat1==sparseLat2,std::string("Both sorts create identical results ")+typeid(data_type).name());

}

BOOST_AUTO_TEST_CASE( SparseLatticeFunctionTests )
{

//do_work<double>(4,4,4);
    do_work<double>(100,100,1000);
    //do_work<double>(2000,2000,10);
//    do_work<double>(100,100,100);
//    do_work<float>(100,100,100);
//    do_work<int>(100,100,100);
//    do_work<long>(100,100,100);



}
