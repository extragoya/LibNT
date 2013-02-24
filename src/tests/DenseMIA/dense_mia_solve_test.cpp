#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE DenseMIASolveTests


#include "DenseMIA.h"
#include "Index.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif


//template<class _data_type>
//void solve_work(size_t dim1){
//
//    LibMIA::MIAINDEX i;
//    LibMIA::MIAINDEX j;
//    LibMIA::MIAINDEX k;
//    LibMIA::MIAINDEX l;
//    LibMIA::MIAINDEX m;
//    LibMIA::MIAINDEX n;
//
//    LibMIA::DenseMIA<_data_type,4> a(dim1,dim1,dim1,dim1);
//    LibMIA::DenseMIA<_data_type,4> b(dim1,dim1,dim1,dim1);
//    LibMIA::DenseMIA<_data_type,4> c(dim1,dim1,dim1,dim1);
//    LibMIA::DenseMIA<_data_type,4> d(dim1,dim1,dim1,dim1);
//
//    a.randu(-30,30);
//    b.randu(-30,30);
//    c(i,j,m,n)=a(i,j,k,l)|b(k,l,m,n);
//    d(i,j,m,n)=a(i,j,k,l)*c(k,l,m,n);
//    BOOST_CHECK_MESSAGE(d.fuzzy_equals(b,test_precision<_data_type>()),std::string("Inner/Outer Product Solve 1 for ")+typeid(_data_type).name() );
//
//
//
//}

BOOST_AUTO_TEST_CASE( DenseMIASolveTests )
{

    std::cout << "AT least we started " << std::endl;
//    solve_work<double>(10);
//    solve_work<float>(10);






}
