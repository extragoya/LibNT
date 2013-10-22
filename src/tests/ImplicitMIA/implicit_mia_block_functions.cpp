#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE ImplicitMIAViewTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "LibMIAUtil.h"
#include "DenseMIA.h"
#include "LibMIARanges.h"
using namespace LibMIA;
template<class _data_type>
void do_work(size_t dim1,size_t dim2){

    DenseMIA<_data_type,2> a(dim1,dim2);
    DenseMIA<_data_type,3> b(dim1,dim2,dim1);
    size_t _start=2,_end=2;
    a.randu(-10,10);
    auto implicit_a=a.view(Range<long long>(_start,dim1-_end),Range<long long>(_start,dim2-_end));
    bool passed=true;


    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            passed=passed&&(implicit_a.at(i,j)==a.at(_start+i,_start+j));
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Value test 1 for view operation for ")+typeid(_data_type).name());
    passed=true;
    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            passed=passed&&(&(implicit_a.at(i,j))==&(a.at(_start+i,_start+j)));
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Reference test 1 for view operation for ")+typeid(_data_type).name());
    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            implicit_a.at(i,j)=2;
        }
    }
    passed=true;
    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            passed=passed&&(a.at(_start+i,_start+j)==2);
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Assign to view test 1 for ")+typeid(_data_type).name());

    //now try with a step size of 2
    _start=1;_end=1;
    size_t _step=2;
    a.randu(-10,10);

    implicit_a=a.view(Range<long long>(_start,dim1-_end,_step),Range<long long>(_start,dim2-_end,_step));

    passed=true;


    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            passed=passed&&(implicit_a.at(i,j)==a.at(_start+i*_step,_start+j*_step));
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Value test 2 for view operation for ")+typeid(_data_type).name());
    passed=true;
    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            passed=passed&&(&(implicit_a.at(i,j))==&(a.at(_start+i*_step,_start+j*_step)));
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Reference test 2 for view operation for ")+typeid(_data_type).name());
    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            implicit_a.at(i,j)=2;
        }
    }
    passed=true;
    for(size_t i=0;i<implicit_a.dim(0);++i){

        for(size_t j=0;j<implicit_a.dim(1);++j){

            passed=passed&&(a.at(_start+i*_step,_start+j*_step)==2);
        }
    }

    BOOST_CHECK_MESSAGE(passed,std::string("Assign to view test 2 for ")+typeid(_data_type).name());

    //now let's try where the we get some 'slices' of the MIA, creating a lower-order MIA

     a.randu(-10,10);

    size_t const_index=3;
    ImplicitMIA<_data_type,1,true> first_order_a=a.view(Range<long long>(_start,dim1-_end,_step),const_index);
    passed=true;


    for(size_t i=0;i<first_order_a.dim(0);++i){
            passed=passed&&(first_order_a.at(i)==a.at(_start+i*_step,const_index));

    }
    BOOST_CHECK_MESSAGE(passed,std::string("Value test 1 for slice for ")+typeid(_data_type).name());

    passed=true;
    for(size_t i=0;i<first_order_a.dim(0);++i){
            passed=passed&&(&(first_order_a.at(i))==&(a.at(_start+i*_step,const_index)));

    }
    BOOST_CHECK_MESSAGE(passed,std::string("Ref test 1 for slice for ")+typeid(_data_type).name());

    for(size_t i=0;i<first_order_a.dim(0);++i){



        first_order_a.at(i)=2;

    }



    passed=true;
    for(size_t i=0;i<first_order_a.dim(0);++i){


        passed=passed&&(a.at(_start+i*_step,const_index)==2);

    }


    BOOST_CHECK_MESSAGE(passed,std::string("Assign to view test 1 for slice for ")+typeid(_data_type).name());


    //Now let's try with other dimension set to constant
    a.randu(-10,10);

    const_index=4;
    first_order_a=a.view(const_index,Range<long long>(_start,dim1-_end,_step));
    passed=true;

    for(size_t i=0;i<first_order_a.dim(0);++i){
            passed=passed&&(first_order_a.at(i)==a.at(const_index,_start+i*_step));

    }
    BOOST_CHECK_MESSAGE(passed,std::string("Value test 2 for slice for ")+typeid(_data_type).name());
    passed=true;
    for(size_t i=0;i<first_order_a.dim(0);++i){
            passed=passed&&(&(first_order_a.at(i))==&(a.at(const_index,_start+i*_step)));

    }
    BOOST_CHECK_MESSAGE(passed,std::string("Ref test 2 for slice for ")+typeid(_data_type).name());

    for(size_t i=0;i<first_order_a.dim(0);++i){



        first_order_a.at(i)=2;

    }
    passed=true;
    for(size_t i=0;i<first_order_a.dim(0);++i){


        passed=passed&&(a.at(const_index,_start+i*_step)==2);

    }

    BOOST_CHECK_MESSAGE(passed,std::string("Assign to view test 2 for slice for ")+typeid(_data_type).name());

    //Now we try with a thiro-order array, using a second-order view
    b.randu(-10,10);
    size_t _step2=4;
    auto implicit_b=b.view(Range<long long>(_start,dim1-_end,_step),const_index,Range<long long>(0,dim1,_step2));

    passed=true;


    for(size_t i=0;i<implicit_b.dim(0);++i){

        for(size_t j=0;j<implicit_b.dim(1);++j){

            passed=passed&&(implicit_b.at(i,j)==b.at(_start+i*_step,const_index,0+j*_step2));
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Value test 3 for slice test for ")+typeid(_data_type).name());
    passed=true;
    for(size_t i=0;i<implicit_b.dim(0);++i){

        for(size_t j=0;j<implicit_b.dim(1);++j){

            passed=passed&&(&(implicit_b.at(i,j))==&(b.at(_start+i*_step,const_index,0+j*_step2)));
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Reference test 3 for view operation for ")+typeid(_data_type).name());
    for(size_t i=0;i<implicit_b.dim(0);++i){

        for(size_t j=0;j<implicit_b.dim(1);++j){

            implicit_b.at(i,j)=2;
        }
    }
    passed=true;
    for(size_t i=0;i<implicit_b.dim(0);++i){

        for(size_t j=0;j<implicit_b.dim(1);++j){

            passed=passed&&(b.at(_start+i*_step,const_index,0+j*_step2)==2);
        }
    }
    BOOST_CHECK_MESSAGE(passed,std::string("Assign to view test 3 for ")+typeid(_data_type).name());



}

BOOST_AUTO_TEST_CASE( ImplicitMIAViewTests )
{



    do_work<double>(8,12);
//    do_work<float>(8,12);
//    do_work<int>(8,12);
//    do_work<long long int>(8,12);




}
