#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE SparseMIAAddSubtractTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#include <boost/timer/timer.hpp>
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


    LibMIA::DenseMIA<_data_type,4> dense_a(dim1,dim1,dim2,dim2);
    LibMIA::DenseMIA<_data_type,4> dense_b(dim2,dim1,dim2,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_c(dim2,dim1,dim2,dim1);

    dense_a.ones();
    dense_b.ones();

    LibMIA::SparseMIA<_data_type,4> a(dense_a);
    LibMIA::SparseMIA<_data_type,4> b(dense_b);
    LibMIA::SparseMIA<_data_type,4> c(dim2,dim1,dim2,dim1);
    //a.print();
    //b.print();
    //boost::timer::cpu_timer scan_t,total_t;
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    //scan_t.stop();
    //c.print();
    //std::cout << "Scan add " << boost::timer::format(scan_t.elapsed()) << std::endl;
    //boost::timer::cpu_timer dense_t;
    dense_c(i,j,k,l)=dense_b(i,j,k,l)+dense_a(j,l,i,k);
    //std::cout << "Dense Scan add " << boost::timer::format(dense_t.elapsed()) << std::endl;
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 1 for ")+typeid(_data_type).name());
    return;

    //test when a is not sorted (uses a different merge algorithm)
    std::array<size_t,4> new_sort_order{{2,0,3,1}};
    a.change_sort_order(new_sort_order);
    //boost::timer::cpu_timer sort_t;
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    //sort_t.stop();
    //std::cout << "Sort add " << boost::timer::format(sort_t.elapsed()) << std::endl;
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Sort) 1 for ")+typeid(_data_type).name());


    a.reset_sort_order();
    a.sort();
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 1 for ")+typeid(_data_type).name());

    a.change_sort_order(new_sort_order);
    dense_b.ones();
    b=dense_b;
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Sort) 1 for ")+typeid(_data_type).name());

    dense_b.init(3);
    b=dense_b;
    a.reset_sort_order();
    a.sort();
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)-dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 1 for ")+typeid(_data_type).name());

    a.change_sort_order(new_sort_order);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Sort) 1 for ")+typeid(_data_type).name());

    a.reset_sort_order();
    a.sort();
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 1 for ")+typeid(_data_type).name());

    a.change_sort_order(new_sort_order);
    dense_b.init(3);
    b=dense_b;
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Sort) 1 for ")+typeid(_data_type).name());

    //****Now try with lots of zeros in our MIAs
    LibMIA::DenseMIA<_data_type,4> dense_a_orig(dim1,dim1,dim2,dim2);
    LibMIA::DenseMIA<_data_type,4> dense_b_orig(dim2,dim1,dim2,dim1);
    dense_a_orig.randu(0,20);
    dense_b_orig.randu(0,20);
    for(auto it=dense_a_orig.data_begin();it<dense_a_orig.data_end();++it)
        if(*it<15)
            *it=0;
    for(auto it=dense_b_orig.data_begin();it<dense_b_orig.data_end();++it)
        if(*it<15)
            *it=0;
    a=dense_a_orig;
    b=dense_b_orig;
    dense_a=dense_a_orig;
    dense_b=dense_b_orig;

    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);


    dense_c(i,j,k,l)=dense_b(i,j,k,l)+dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 2 for ")+typeid(_data_type).name());


    //test when a is not sorted (uses a different merge algorithm)
    a.change_sort_order(new_sort_order);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Sort) 2 for ")+typeid(_data_type).name());


    a.reset_sort_order();
    a.sort();
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 2 for ")+typeid(_data_type).name());

    a.change_sort_order(new_sort_order);
    dense_b=dense_b_orig;
    b=dense_b;
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Sort) 2 for ")+typeid(_data_type).name());

    dense_b=dense_b_orig;
    b=dense_b;
    a.reset_sort_order();

    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)-dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 2 for ")+typeid(_data_type).name());

    a.change_sort_order(new_sort_order);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Sort) 2 for ")+typeid(_data_type).name());

    a.reset_sort_order();
    a.sort();
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 2 for ")+typeid(_data_type).name());

    a.change_sort_order(new_sort_order);
    dense_b=dense_b_orig;
    b=dense_b;
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Sort) 2 for ")+typeid(_data_type).name());



}

BOOST_AUTO_TEST_CASE( SparseMIAAddSubtractTests )
{



    do_work<double>(10,80);
//    do_work<double>(10,15);
//    do_work<float>(10,15);
//    do_work<int>(10,15);
//    do_work<long long int>(10,15);

//    do_work<double>(10,8);
//    do_work<float>(10,8);
//    do_work<int>(10,8);
//    do_work<long long int>(10,8);




}
