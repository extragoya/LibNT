#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE MixedMIAMultTests
#include "MIAConfig.h"
#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "SparseMIA.h"
#include "DenseMIA.h"
#include "Index.h"





template<class _data_type>
void mult_work(size_t dim1, size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;

    LibMIA::DenseMIA<_data_type,4> temp_a(dim1,dim2,dim1,dim2);
    LibMIA::DenseMIA<_data_type,4> dense_b(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_c;
    LibMIA::DenseMIA<_data_type,2> dense_c2;
    LibMIA::DenseMIA<_data_type,2> dense_b2(dim2,dim2);
    LibMIA::DenseMIA<_data_type,1> dense_d(dim2);
    LibMIA::DenseMIA<_data_type,2> dense_d2(dim2,dim2);
    LibMIA::DenseMIA<_data_type,4> c_mixed;
    LibMIA::DenseMIA<_data_type,2> c2_mixed;

    LibMIA::SparseMIA<_data_type,4> a(dim1,dim2,dim1,dim2);
    LibMIA::SparseMIA<_data_type,4> b(dim2,dim2,dim1,dim1);

    LibMIA::SparseMIA<_data_type,2> b2(dim2,dim2);
    LibMIA::SparseMIA<_data_type,1> d(dim2);
    LibMIA::SparseMIA<_data_type,2> d2(dim2,dim2);
    LibMIA::SparseMIA<_data_type,4> c;


    temp_a.randu(0,20);
    dense_b.randu(0,20);
    dense_b2.randu(0,20);
    dense_d.randu(0,20);
    dense_d2.randu(0,20);

    for(auto it=temp_a.data_begin();it<temp_a.data_end();++it)
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

    const LibMIA::DenseMIA<_data_type,4> dense_a(temp_a);
    a=dense_a;
    b=dense_b;
    b2=dense_b2;
    d=dense_d;
    d2=dense_d2;

    dense_c(i,k,m,n)=dense_a(i,j,k,l)*dense_b(j,l,m,n);
    c_mixed(i,k,m,n)=a(i,j,k,l)*dense_b(j,l,m,n);
    BOOST_CHECK_MESSAGE(c_mixed.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product 1a for ")+typeid(_data_type).name() );


    c_mixed(i,k,m,n)=dense_a(i,j,k,l)*b(j,l,m,n);
    BOOST_CHECK_MESSAGE(c_mixed.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product 1b for ")+typeid(_data_type).name() );


    dense_c(i,k,m,n)=dense_a(i,l,k,j)*dense_b(l,j,m,n);
    c_mixed(i,k,m,n)=a(i,l,k,j)*dense_b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c_mixed.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product 2a for ")+typeid(_data_type).name());
    c_mixed(i,k,m,n)=dense_a(i,l,k,j)*b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c_mixed.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product 2b for ")+typeid(_data_type).name());

    dense_c(i,k,m,n)=dense_a(i,l,k,j)*dense_b(l,j,m,n);
    c_mixed(i,k,m,n)=a(i,l,k,j)*dense_b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c_mixed.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product 3a for ")+typeid(_data_type).name());
    c_mixed(i,k,m,n)=dense_a(i,l,k,j)*b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c_mixed.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product 3b for ")+typeid(_data_type).name());


    dense_c2(i,j)=dense_a(!i,k,!j,l)*dense_b(k,l,!i,!j);
    c2_mixed(i,j)=a(!i,k,!j,l)*dense_b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2_mixed.fuzzy_equals(dense_c2,test_precision<_data_type>()),std::string("Inner/Element-Wise Product 1a for ")+typeid(_data_type).name());
    c2_mixed(i,j)=dense_a(!i,k,!j,l)*b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2_mixed.fuzzy_equals(dense_c2,test_precision<_data_type>()),std::string("Inner/Element-Wise Product 1b for ")+typeid(_data_type).name());

    dense_c2(i,j)=dense_a(!i,l,!j,k)*dense_b(k,l,!i,!j);
    c2_mixed(i,j)=a(!i,l,!j,k)*dense_b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2_mixed.fuzzy_equals(dense_c2,test_precision<_data_type>()),std::string("Inner/Element-Wise Product 2a for ")+typeid(_data_type).name());
    c2_mixed(i,j)=dense_a(!i,l,!j,k)*b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2_mixed.fuzzy_equals(dense_c2,test_precision<_data_type>()),std::string("Inner/Element-Wise Product 2b for ")+typeid(_data_type).name());

    dense_c2(j,i)=dense_a(!j,l,!i,k)*dense_b(k,l,!j,!i);
    c2_mixed(j,i)=a(!j,l,!i,k)*dense_b(k,l,!j,!i);
    BOOST_CHECK_MESSAGE(c2_mixed.fuzzy_equals(dense_c2,test_precision<_data_type>()),std::string("Inner/Element-Wise Product 3a for ")+typeid(_data_type).name());
    c2_mixed(j,i)=dense_a(!j,l,!i,k)*b(k,l,!j,!i);
    BOOST_CHECK_MESSAGE(c2_mixed.fuzzy_equals(dense_c2,test_precision<_data_type>()),std::string("Inner/Element-Wise Product 3b for ")+typeid(_data_type).name());


    dense_c(i,j,k,l)=dense_b2(i,j)*dense_d2(k,l);
    c(i,j,k,l)=b2(i,j)*dense_d2(k,l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Outer Product 1a for ")+typeid(_data_type).name());
    c(i,j,k,l)=dense_b2(i,j)*d2(k,l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Outer Product 1b for ")+typeid(_data_type).name());


    dense_c(i,k,l,j)=dense_b2(k,j)*dense_d2(l,i);
    c(i,k,l,j)=b2(k,j)*dense_d2(l,i);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Outer Product 2a for ")+typeid(_data_type).name());
    c(i,k,l,j)=dense_b2(k,j)*d2(l,i);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Outer Product 2b for ")+typeid(_data_type).name());

    dense_c(i,k,l,j)=dense_d2(l,i)*dense_b2(k,j);
    c(i,k,l,j)=d2(l,i)*dense_b2(k,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Outer Product 3a for ")+typeid(_data_type).name());
    c(i,k,l,j)=dense_d2(l,i)*b2(k,j);

//    dense_c.print();
//    c.print();
//    c.reset_linIdx_sequence();
//    c.print();
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Outer Product 3b for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=dense_a(i,!j,k,!l)*dense_b2(!j,!l);
    c(i,j,k,l)=a(i,!j,k,!l)*dense_b2(!j,!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Element-Wise Product 1a for ")+typeid(_data_type).name());
    c(i,j,k,l)=dense_a(i,!j,k,!l)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Element-Wise Product 1b for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=dense_a(k,!j,i,!l)*dense_b2(!j,!l);
    c(i,j,k,l)=a(k,!j,i,!l)*dense_b2(!j,!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Element-Wise Product 2a for ")+typeid(_data_type).name());
    c(i,j,k,l)=dense_a(k,!j,i,!l)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Element-Wise Product 2b for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=dense_a(k,!l,i,!j)*dense_b2(!j,!l);
    c(i,j,k,l)=a(k,!l,i,!j)*dense_b2(!j,!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Element-Wise Product 3a for ")+typeid(_data_type).name());
    c(i,j,k,l)=dense_a(k,!l,i,!j)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Outer/Element-Wise Product 3b for ")+typeid(_data_type).name());

    dense_c(i,j,k,l)=~(dense_a(i,!j,k,!!l)*dense_b2(!j,!!l))*dense_d(!l);
    c(i,j,k,l)=~(a(i,!j,k,!!l)*dense_b2(!j,!!l))*d(!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Repeated Element-Wise Product 1a for ")+typeid(_data_type).name());
    c(i,j,k,l)=~(dense_a(i,!j,k,!!l)*b2(!j,!!l))*d(!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Repeated Element-Wise Product 1b for ")+typeid(_data_type).name());
    c(i,j,k,l)=~(a(i,!j,k,!!l)*b2(!j,!!l))*dense_d(!l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Repeated Element-Wise Product 1c for ")+typeid(_data_type).name());


    dense_c(i,k,m,n)=~(dense_a(i,!j,k,!l)*dense_b(!j,!l,m,n))*dense_d2(j,l);

    c(i,k,m,n)=~(a(i,!j,k,!l)*dense_b(!j,!l,m,n))*d2(j,l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 1a for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(dense_a(i,!j,k,!l)*b(!j,!l,m,n))*d2(j,l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 1b for ")+typeid(_data_type).name() );

    c(i,k,m,n)=~(a(i,!j,k,!l)*b(!j,!l,m,n))*dense_d2(j,l);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 1c for ")+typeid(_data_type).name() );

    dense_c(i,k,m,n)=~(dense_a(i,!j,k,!l)*dense_b(!j,!l,m,n))*dense_d2(l,j);
    c(i,k,m,n)=~(a(i,!j,k,!l)*dense_b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 2a for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(dense_a(i,!j,k,!l)*b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 2b for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(a(i,!j,k,!l)*b(!j,!l,m,n))*dense_d2(l,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 2c for ")+typeid(_data_type).name() );

    dense_c(i,k,m,n)=~(dense_a(i,!l,k,!j)*dense_b(!j,!l,m,n))*dense_d2(l,j);
    c(i,k,m,n)=~(a(i,!l,k,!j)*dense_b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 3a for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(dense_a(i,!l,k,!j)*b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 3b for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(a(i,!l,k,!j)*b(!j,!l,m,n))*dense_d2(l,j);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Ternary Inner Product 3c for ")+typeid(_data_type).name() );


}

BOOST_AUTO_TEST_CASE( MixedMIAMultTests )
{

    //mult_work<double>(3,3);
    //mult_work<float>(3,3);
    //mult_work<int>(3,3);
    //mult_work<long>(3,3);
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
