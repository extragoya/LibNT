#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE DenseMIAAddSubtractTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "DenseMIA.h"
#include "Index.h"

template<class _data_type>
void do_work(size_t dim1,size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
//    LibMIA::MIAINDEX m;
//    LibMIA::MIAINDEX n;
//    LibMIA::MIAINDEX o;
//    LibMIA::MIAINDEX p;

    LibMIA::DenseMIA<_data_type,4> a(dim1,dim1,dim2,dim2);
    LibMIA::DenseMIA<_data_type,4> b(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> c(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> c2(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> b2(dim1,dim1,dim2,dim2);


    a.ones();
    b.ones();


    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    c2.init(2);
    BOOST_CHECK_MESSAGE(c==c2,std::string("Non-destructive Add 1 for ")+typeid(_data_type).name());


    b(i,j,k,l)+=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==c2,std::string("Destructive Add 1 for ")+typeid(_data_type).name());

    b.init(3);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==c2,std::string("Non-destructive Subtract 1 for ")+typeid(_data_type).name());


    b(i,j,k,l)-=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==c2,std::string("Destructive Subtract 1 for ")+typeid(_data_type).name());



}

BOOST_AUTO_TEST_CASE( DenseMIAAddSubtractTests )
{



    do_work<double>(8,10);
    do_work<float>(8,10);
    do_work<int>(8,10);
    do_work<long long int>(8,10);




}
