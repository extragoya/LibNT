#define BOOST_TEST_MODULE SparseLatticeMultTestsSpecial




#include "SparseLattice.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

//!Relies on normal SparseLattice*SparseLattice operator functioning correctly, i.e., that the SparseLatticeMultTests passed
template<typename data_type>
void multwork(size_t m1, size_t n1, size_t n2, size_t p,double hypersparsity){




    typedef LibMIA::SparseLattice<data_type> sparseType;
    sparseType A,B;
    A=sparseType(m1,n1,p);
    A.resize(std::ceil(p*m1*hypersparsity));
    A.randu(1,10);
    A.rand_indices();
    A.collect_duplicates();

    B=sparseType(n1,n2,p);
    B.resize(std::ceil(p*n2*hypersparsity));
    B.randu(1,10);
    B.rand_indices();
    B.collect_duplicates();

    sparseType C,C_test;

    C=A*B;

    C_test=A.template csc_times<false>(B);
    BOOST_CHECK_MESSAGE(C.fuzzy_equals(C_test,test_precision<data_type>()),std::string("CSC Mult Test Accum for ")+typeid(data_type).name());
	

    C_test=A.template csc_no_accum<false>(B);
    BOOST_CHECK_MESSAGE(C.fuzzy_equals(C_test,test_precision<data_type>()),std::string("CSC Mult Test No Accum for ")+typeid(data_type).name());
	
    C_test=A.template csc_times<true>(B);
    BOOST_CHECK_MESSAGE(C.fuzzy_equals(C_test,test_precision<data_type>()),std::string("DCSC Mult Test Accum for ")+typeid(data_type).name());

    C_test=A.template csc_no_accum<true>(B);
    BOOST_CHECK_MESSAGE(C.fuzzy_equals(C_test,test_precision<data_type>()),std::string("DCSC Mult Test No Accum for ")+typeid(data_type).name());

    C_test=A.outer_times(B);
    BOOST_CHECK_MESSAGE(C.fuzzy_equals(C_test,test_precision<data_type>()),std::string("Outer Mult Test for ")+typeid(data_type).name());




}

BOOST_AUTO_TEST_CASE( SparseLatticeMultTestsSpecial )
{


    //multwork<double>(5,5,5,5,1);
    multwork<double>(20,20,20,20,1);
    multwork<float>(20,20,20,20,1);
    multwork<int>(20,20,20,20,1);
    multwork<long>(20,20,20,20,1);

    multwork<double>(20,20,20,20,.5);
    multwork<float>(20,20,20,20,.5);
    multwork<int>(20,20,20,20,.5);
    multwork<long>(20,20,20,20,.5);

    multwork<double>(40,20,20,20,1);
    multwork<float>(40,20,20,20,1);
    multwork<int>(40,20,20,20,1);
    multwork<long>(40,20,20,20,1);

    multwork<double>(20,20,40,20,.5);
    multwork<float>(20,20,40,20,.5);
    multwork<int>(20,20,20,40,.5);
    multwork<long>(20,20,20,40,.5);



}
