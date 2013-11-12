#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE SparseMIAMultTests
#include "MIAConfig.h"
#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "SparseMIA.h"
#include "DenseMIA.h"
#include "Index.h"
#include "LibMIAUtil.h"
#include "FunctionUtil.h"



template<class _data_type>
void mult_work(size_t dim1, size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;

    LibMIA::DenseMIA<_data_type,4> dense_a(dim1,dim2,dim1,dim2);
    LibMIA::DenseMIA<_data_type,4> dense_b(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_c;
    LibMIA::DenseMIA<_data_type,2> dense_c2;
    LibMIA::DenseMIA<_data_type,3> dense_c3;
    LibMIA::DenseMIA<_data_type,2> dense_b2(dim2,dim2);
    LibMIA::DenseMIA<_data_type,1> dense_d(dim2);
    LibMIA::DenseMIA<_data_type,2> dense_d2(dim2,dim2);

    LibMIA::SparseMIA<_data_type,4> a(dim1,dim2,dim1,dim2);
    LibMIA::SparseMIA<_data_type,4> b(dim2,dim2,dim1,dim1);
    LibMIA::SparseMIA<_data_type,4> c;
    LibMIA::SparseMIA<_data_type,2> c2;
    LibMIA::SparseMIA<_data_type,3> c3;
    LibMIA::SparseMIA<_data_type,2> b2(dim2,dim2);
    LibMIA::SparseMIA<_data_type,1> d(dim2);
    LibMIA::SparseMIA<_data_type,2> d2(dim2,dim2);


    dense_a.randu(0,20);
    dense_b.randu(0,20);
    dense_b2.randu(0,20);
    dense_d.randu(0,20);
    dense_d2.randu(0,20);

    for(auto it=dense_a.data_begin();it<dense_a.data_end();++it)
        if(*it<15)
            *it=0;
    for(auto it=dense_b.data_begin();it<dense_b.data_end();++it)
        if(*it<15)
            *it=0;
    for(auto it=dense_b2.data_begin();it<dense_b2.data_end();++it)
        if(*it<15)
            *it=0;
    for(auto it=dense_d.data_begin();it<dense_d.data_end();++it)
        if(*it<15)
            *it=0;
    for(auto it=dense_d2.data_begin();it<dense_d2.data_end();++it)
        if(*it<15)
            *it=0;

    a=dense_a;
    b=dense_b;
    b2=dense_b2;
    d=dense_d;
    d2=dense_d2;

    dense_c(i,k,m,n)=dense_a(i,j,k,l)*dense_b(j,l,m,n);
    c(i,k,m,n)=a(i,j,k,l)*b(j,l,m,n);

    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Inner/Outer Product 1 for ")+typeid(_data_type).name() );


    dense_c(l,k,m,n)=dense_a(i,j,m,n)*dense_a(k,l,i,j);
    c(l,k,m,n)=a(i,j,m,n)*a(k,l,i,j);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Inner/Outer Product 1 with self-multiplication for ")+typeid(_data_type).name() ); //test the self-multiplication

    dense_c(i,k,m,n)=dense_a(i,l,k,j)*dense_b(l,j,m,n);
    c(i,k,m,n)=a(i,l,k,j)*b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Inner/Outer Product 2 for ")+typeid(_data_type).name());

    dense_c(i,k,m,n)=dense_a(i,l,k,j)*dense_b(l,j,m,n);
    c(i,k,m,n)=a(i,l,k,j)*b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Inner/Outer Product 3 for ")+typeid(_data_type).name());


    dense_c2(i,j)=dense_a(!i,k,!j,l)*dense_b(k,l,!i,!j);
    c2(i,j)=a(!i,k,!j,l)*b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2==dense_c2,std::string("Inner/Element-Wise Product 1 for ")+typeid(_data_type).name());

    dense_c2(i,j)=dense_a(!i,l,!j,k)*dense_b(k,l,!i,!j);
    c2(i,j)=a(!i,l,!j,k)*b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2==dense_c2,std::string("Inner/Element-Wise Product 2 for ")+typeid(_data_type).name());

    dense_c2(j,i)=dense_a(!j,l,!i,k)*dense_b(k,l,!j,!i);
    c2(j,i)=a(!j,l,!i,k)*b(k,l,!j,!i);
    BOOST_CHECK_MESSAGE(c2==dense_c2,std::string("Inner/Element-Wise Product 3 for ")+typeid(_data_type).name());


    dense_c(i,j,k,l)=dense_b2(i,j)*dense_d2(k,l);
    c(i,j,k,l)=b2(i,j)*d2(k,l);
    //dense_c.print();
    //c.print();
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Outer/Outer Product 1 for ")+typeid(_data_type).name());

    dense_c(i,k,l,j)=dense_b2(k,j)*dense_d2(l,i);
    c(i,k,l,j)=b2(k,j)*d2(l,i);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Outer/Outer Product 2 for ")+typeid(_data_type).name());

    dense_c(i,k,l,j)=dense_d2(l,i)*dense_b2(k,j);
    c(i,k,l,j)=d2(l,i)*b2(k,j);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Outer/Outer Product 3 for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=dense_a(i,!j,k,!l)*dense_b2(!j,!l);
    c(i,j,k,l)=a(i,!j,k,!l)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Outer/Element-Wise Product 1 for ")+typeid(_data_type).name());

    dense_c3(i,j,l)=dense_b2(i,!j)*dense_b2(!j,l);
    c3(i,j,l)=b2(i,!j)*b2(!j,l);
    BOOST_CHECK_MESSAGE(c3==dense_c3,std::string("Outer/Element-Wise Product 1 with self-multiplication for ")+typeid(_data_type).name()); //test the self-multiplication

    dense_c(i,j,k,l)=dense_a(k,!j,i,!l)*dense_b2(!j,!l);
    c(i,j,k,l)=a(k,!j,i,!l)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Outer/Element-Wise Product 2 for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=dense_a(k,!l,i,!j)*dense_b2(!j,!l);
    c(i,j,k,l)=a(k,!l,i,!j)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Outer/Element-Wise Product 3 for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=~(dense_a(i,!j,k,!!l)*dense_b2(!j,!!l))*dense_d(!l);
    c(i,j,k,l)=~(a(i,!j,k,!!l)*b2(!j,!!l))*d(!l);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Repeated Element-Wise Product 1 for ")+typeid(_data_type).name());


    dense_c(i,k,m,n)=~(dense_a(i,!j,k,!l)*dense_b(!j,!l,m,n))*dense_d2(j,l);
    c(i,k,m,n)=~(a(i,!j,k,!l)*b(!j,!l,m,n))*d2(j,l);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Ternary Inner Product 1 for ")+typeid(_data_type).name() );

    dense_c(i,k,m,n)=~(dense_a(i,!j,k,!l)*dense_b(!j,!l,m,n))*dense_d2(l,j);
    c(i,k,m,n)=~(a(i,!j,k,!l)*b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Ternary Inner Product 2 for ")+typeid(_data_type).name() );

    dense_c(i,k,m,n)=~(dense_a(i,!l,k,!j)*dense_b(!j,!l,m,n))*dense_d2(l,j);
    c(i,k,m,n)=~(a(i,!l,k,!j)*b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Ternary Inner Product 3 for ")+typeid(_data_type).name() );

    _data_type dense_data=dense_a(i,l,k,j)*dense_b(i,l,j,k);
    _data_type sparse_data=a(i,l,k,j)*b(i,l,j,k);
    if(dense_data!=sparse_data)
        std::cout << "dense " << dense_data << " sparse_data " << sparse_data << std::endl;
    BOOST_CHECK_MESSAGE(LibMIA::isEqualFuzzy(dense_data,sparse_data,test_precision<_data_type>()),std::string("Complete Inner Product test 1 for ")+typeid(_data_type).name() );

    dense_data=dense_a(i,l,k,j)*dense_a(i,j,k,l);
    sparse_data=a(i,l,k,j)*a(i,j,k,l);
    BOOST_CHECK_MESSAGE(LibMIA::isEqualFuzzy(dense_data,sparse_data,test_precision<_data_type>()),std::string("Complete Inner Product test 1 with self-assignment for ")+typeid(_data_type).name() );

}

BOOST_AUTO_TEST_CASE( SparseMIAMultTests )
{

    //mult_work<double>(3,4);
//    mult_work<float>(3,3);
//    mult_work<int>(3,3);
//    mult_work<long>(3,3);
//
//


    mult_work<double>(8,5);
    mult_work<float>(8,5);
    mult_work<int>(8,5);
    mult_work<long>(8,5);



    mult_work<double>(5,8);
    mult_work<float>(5,8);
    mult_work<int>(5,8);
    mult_work<long>(5,8);


}
