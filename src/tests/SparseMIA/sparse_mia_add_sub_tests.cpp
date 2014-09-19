#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE SparseMIAAddSubtractTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>

#else
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#endif



#include "DenseMIA.h"
#include "SparseMIA.h"
#include "Index.h"
#include "LibMIAUtil.h"
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

    std::array<size_t,4> new_linIdxSequence{{2,0,3,1}};
    std::array<size_t,4> new_linIdxSequence2{{1,3,0,2}};
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);

    //scan_t.stop();
    //c.print();
    //std::cout << "Scan add " << boost::timer::format(scan_t.elapsed()) << std::endl;
    //boost::timer::cpu_timer dense_t;
    dense_c(i,j,k,l)=dense_b(i,j,k,l)+dense_a(j,l,i,k);
    //std::cout << "Dense Scan add " << boost::timer::format(dense_t.elapsed()) << std::endl;

    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 1a for ")+typeid(_data_type).name());

    //now check when b is not sorted
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 1b for ")+typeid(_data_type).name());

    //now check when a and b have non default sort orders
    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 1c for ")+typeid(_data_type).name());


    c(i,j,k,l)=b(i,j,k,l)+b(k,l,i,j);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)+dense_b(k,l,i,j);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Self addition check 1 for ")+typeid(_data_type).name());


    //test when a is not sorted (uses a different merge algorithm)
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b.change_linIdx_sequence(new_linIdxSequence);
    //boost::timer::cpu_timer sort_t;
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);

    //sort_t.stop();
    //std::cout << "Sort add " << boost::timer::format(sort_t.elapsed()) << std::endl;
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Sort) 1a for ")+typeid(_data_type).name());

    //test when a is not sorted (uses a different merge algorithm) and both have different sort orders
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b.change_linIdx_sequence(new_linIdxSequence2);
    //boost::timer::cpu_timer sort_t;
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);

    //sort_t.stop();
    //std::cout << "Sort add " << boost::timer::format(sort_t.elapsed()) << std::endl;
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Sort) 1b for ")+typeid(_data_type).name());



    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 1a for ")+typeid(_data_type).name());

    dense_b.ones();
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 1b for ")+typeid(_data_type).name());

    dense_b.ones();
    b=dense_b;
    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 1c for ")+typeid(_data_type).name());


    dense_b.ones();
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Sort) 1a for ")+typeid(_data_type).name());

    dense_b.ones();
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b.change_linIdx_sequence(new_linIdxSequence2);
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Sort) 1b for ")+typeid(_data_type).name());

    dense_b.fill(3);
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)-dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 1a for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 1b for ")+typeid(_data_type).name());

    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 1c for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Sort) 1a for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    a.change_linIdx_sequence(new_linIdxSequence2);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Sort) 1b for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 1a for ")+typeid(_data_type).name());

    dense_b.fill(3);
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 1b for ")+typeid(_data_type).name());

    dense_b.fill(3);
    b=dense_b;
    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 1c for ")+typeid(_data_type).name());

    dense_b.fill(3);
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Sort) 1a for ")+typeid(_data_type).name());

    dense_b.fill(3);
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence2);
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Sort) 1b for ")+typeid(_data_type).name());

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

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)+dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 2a for ")+typeid(_data_type).name());



    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 2b for ")+typeid(_data_type).name());

    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Scan) 2c for ")+typeid(_data_type).name());

    //test when a is not sorted (uses a different merge algorithm)
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Sort) 2a for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Add (Sort) 2b for ")+typeid(_data_type).name());


    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)+a(j,l,i,k)+b(k,l,i,j)+a(l,j,i,k);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)+dense_a(j,l,i,k)+dense_b(k,l,i,j)+dense_a(l,j,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Chained Non-destructive Add (Sort) 2a for ")+typeid(_data_type).name());


    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    dense_b(i,j,k,l)+=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 2a for ")+typeid(_data_type).name());

    b=dense_b_orig;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 2b for ")+typeid(_data_type).name());

    b=dense_b_orig;
    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    b(i,j,k,l)+=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Scan) 2c for ")+typeid(_data_type).name());


    b=dense_b_orig;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Sort) 2a for ")+typeid(_data_type).name());

    b=dense_b_orig;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)+=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Add (Sort) 2b for ")+typeid(_data_type).name());

    dense_b=dense_b_orig;
    b=dense_b;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    dense_c(i,j,k,l)=dense_b(i,j,k,l)-dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 2a for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 2b for ")+typeid(_data_type).name());

    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Scan) 2c for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Sort) 2a for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    a.change_linIdx_sequence(new_linIdxSequence);
    c(i,j,k,l)=b(i,j,k,l)-a(j,l,i,k);
    BOOST_CHECK_MESSAGE(c==dense_c,std::string("Non-destructive Subtract (Sort) 2b for ")+typeid(_data_type).name());

    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    dense_b(i,j,k,l)-=dense_a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 2a for ")+typeid(_data_type).name());



    b=dense_b_orig;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 2b for ")+typeid(_data_type).name());

    b=dense_b_orig;
    a.change_linIdx_sequence(new_linIdxSequence);
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    b(i,j,k,l)-=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Scan) 2c for ")+typeid(_data_type).name());


    b=dense_b_orig;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence);
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Sort) 2a for ")+typeid(_data_type).name());

    b=dense_b_orig;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();
    b.change_linIdx_sequence(new_linIdxSequence2);
    a.change_linIdx_sequence(new_linIdxSequence);
    b(i,j,k,l)-=a(j,l,i,k);
    BOOST_CHECK_MESSAGE(b==dense_b,std::string("Destructive Subtract (Sort) 2b for ")+typeid(_data_type).name());



}

BOOST_AUTO_TEST_CASE( SparseMIAAddSubtractTests )
{



    //do_work<double>(2,3);
    do_work<double>(5,8);
    do_work<float>(5,8);
    do_work<int>(5,8);
    do_work<long long int>(5,8);

//    do_work<double>(10,8);
//    do_work<float>(10,8);
//    do_work<int>(10,8);
//    do_work<long long int>(10,8);




}
