#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE SparseMIAAddSubtractTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"
#include "SparseMIA.h"
#include "Index.h"
#include "Util.h"
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



    LibMIA::DenseMIA<_data_type,4> dense_a(dim1,dim1,dim2,dim2);
    LibMIA::DenseMIA<_data_type,4> dense_b(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_c(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_c2(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_b2(dim1,dim1,dim2,dim2);

    dense_a.ones();
    dense_b.ones();


    //c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    dense_c2.init(2);
    //BOOST_CHECK_MESSAGE(c==c2,std::string("Non-destructive Add 1 for ")+typeid(_data_type).name());
    LibMIA::SparseMIA<_data_type,4> a(dense_a);

    LibMIA::SparseMIA<_data_type,4> b(dense_b);
    std::cout << "**********A*************" <<std::endl;
    a.print();
    std::cout << "**********B*************" <<std::endl;
    a.print();
    b(i,j,k,l)+=a(j,l,i,k);
    std::cout << "**********B RESULT*************" <<std::endl;
    b.print();
    BOOST_CHECK_MESSAGE(b==dense_c2,std::string("Destructive Add 1 for ")+typeid(_data_type).name());

//    b.init(3);
//    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
//    BOOST_CHECK_MESSAGE(c==c2,std::string("Non-destructive Subtract 1 for ")+typeid(_data_type).name());
//
//
//    b(i,j,k,l)-=a(j,l,i,k);
//    BOOST_CHECK_MESSAGE(b==c2,std::string("Destructive Subtract 1 for ")+typeid(_data_type).name());
//
//
//    a.zeros();
//    b.zeros();
//    a.at(dim1-1,dim1-1,dim2-1,dim2-1)=1;
//    b.at(dim2-1,dim1-1,dim2-1,dim1-1)=1;
//
//
//    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
//    c2.zeros();
//    c2.at(dim2-1,dim1-1,dim2-1,dim1-1)=2;
//    BOOST_CHECK_MESSAGE(c==c2,std::string("Non-destructive Add 2 for ")+typeid(_data_type).name());
//
//
//    b(i,j,k,l)+=a(j,l,i,k);
//    BOOST_CHECK_MESSAGE(b==c2,std::string("Destructive Add 2 for ")+typeid(_data_type).name());
//
//    a.at(dim1-1,dim1-1,dim2-1,dim2-1)=3;
//    b.at(dim2-1,dim1-1,dim2-1,dim1-1)=5;
//    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
//    BOOST_CHECK_MESSAGE(c==c2,std::string("Non-destructive Subtract 2 for ")+typeid(_data_type).name());
//
//
//    b(i,j,k,l)-=a(j,l,i,k);
//    BOOST_CHECK_MESSAGE(b==c2,std::string("Destructive Subtract 2 for ")+typeid(_data_type).name());



}

BOOST_AUTO_TEST_CASE( SparseMIAAddSubtractTests )
{

    do_work<double>(2,3);


    //do_work<double>(8,10);
//    do_work<float>(8,10);
//    do_work<int>(8,10);
//    do_work<long long int>(8,10);
//
//    do_work<double>(10,8);
//    do_work<float>(10,8);
//    do_work<int>(10,8);
//    do_work<long long int>(10,8);




}
