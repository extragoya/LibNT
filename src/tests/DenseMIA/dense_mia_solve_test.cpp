#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE DenseMIASolveTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "DenseMIA.h"
#include "Index.h"

template<class _data_type>
void solve_work(size_t dim1,size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;
    LibMIA::MIAINDEX o;
    LibMIA::MIAINDEX p;

    LibMIA::DenseMIA<_data_type,4> a(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> a2(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> b(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> b2(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> c(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> d(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> d2(dim1,dim1,dim1,dim1);

    a.randu(-30,30); //we can invert these, b/c random matrices are non-singular with probability one (almost surely)
    b.randu(-30,30);
    a2.randu(-30,30);
    b2.randu(-30,30);

    c(i,j,m,n)=a(i,j,k,l)|b(k,l,m,n);
    d(k,l,m,n)=a(i,j,k,l)*c(i,j,m,n);
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(b,test_precision<_data_type>()),std::string("Inner/Outer Product Inverse 1 for ")+typeid(_data_type).name() );

    c(i,j,m,n)=a(i,k,j,l)|b(k,l,m,n);
    d(k,l,m,n)=a(i,k,j,l)*c(i,j,m,n);
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(b,test_precision<_data_type>()),std::string("Inner/Outer Product Inverse 2 for ")+typeid(_data_type).name() );

    c(i,j,l,m)=a(i,!j,k,!l)|b(k,!j,!l,m);
    d(k,j,l,m)=a(i,!j,k,!l)*c(i,!j,!l,m);
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(b,test_precision<_data_type>()),std::string("Inner/Outer/Inter Product Inverse 1 for ")+typeid(_data_type).name() );

    c(i,j,l,m)=a(!i,!j,k,l)|b(m,!j,k,!i);
    d(m,j,k,i)=a(!i,!j,k,l)*c(!i,!j,l,m);
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(b,test_precision<_data_type>()),std::string("Inner/Outer/Inter Product Inverse 2 for ")+typeid(_data_type).name() );

    c(i,j,m,n)=a2(k,l,i,j)|b2(k,l,m,n);
    //test with normal equations
    d(o,p,m,n)=a2(k,l,o,p)*a2(k,l,i,j)*c(i,j,m,n);
    d2(i,j,m,n)=a2(k,l,i,j)*b2(k,l,m,n);
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(d2,test_precision<_data_type>()),std::string("Inner/Outer Product Least Squares 1 for ")+typeid(_data_type).name() );

}

BOOST_AUTO_TEST_CASE( DenseMIASolveTests )
{



    solve_work<double>(8,10);
    solve_work<float>(8,10);






}
