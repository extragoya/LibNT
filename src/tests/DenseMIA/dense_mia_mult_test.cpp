#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE DenseMIAMultTests


#include "DenseMIA.h"
#include "Index.h"
#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif


template<class _data_type>
void mult_work(size_t dim1, size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;

    LibMIA::DenseMIA<_data_type,4> a(dim1,dim2,dim1,dim2);
    LibMIA::DenseMIA<_data_type,4> b(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> c;
    LibMIA::DenseMIA<_data_type,2> c2;
    LibMIA::DenseMIA<_data_type,4> c_result(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,2> c_result2(dim1,dim1);

    b.ones();

    //set each the values of a to increment along inner product indices
    size_t val,entries=dim2*dim2;
    for(size_t _i=0;_i<dim1;++_i){
        for(size_t _k=0;_k<dim1;++_k){
            val=1;
            for(size_t _j=0;_j<dim2;++_j){
                for(size_t _l=0;_l<dim2;++_l){
                    a.at(_i,_j,_k,_l)=val++;
                }
            }
        }
    }

    for(size_t _i=0;_i<dim1;++_i){
        for(size_t _j=0;_j<dim1;++_j){

            for(size_t _k=0;_k<dim1;++_k){
                for(size_t _l=0;_l<dim1;++_l){
                    c_result.at(_i,_j,_k,_l)=entries;
                }
            }
        }
    }

    c(i,k,m,n)=a(i,j,k,l)*b(j,l,m,n);
    BOOST_CHECK_MESSAGE(c==c_result,std::string("Inner/Outer Product 1 for ")+typeid(_data_type).name() );

    c(i,k,m,n)=a(i,l,k,j)*b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c==c_result,std::string("Inner/Outer Product 2 for ")+typeid(_data_type).name());

    c(i,k,m,n)=a(i,j,k,l)*b(l,j,m,n);
    BOOST_CHECK_MESSAGE(c==c_result,std::string("Inner/Outer Product 3 for ")+typeid(_data_type).name());

//    //test inner and element-wise product. Have b be assigned a scalar value that increases while traversing i
//    val=0;
//    for(size_t _i=0;_i<dim1;++_i){
//        val++;
//        for(size_t _j=0;_j<dim1;++_j){
//            c_result2(i,j)=val*entries;
//            for(size_t _k=0;_k<dim2;++_k){
//                for(size_t _l=0;_l<dim2;++_l){
//                    b.at(_k,_l,_i,_j)=val;
//                }
//            }
//        }
//    }
//
//    c2(i,j)=a(!i,k,!j,l)*b(k,l,!i,!j);
//    BOOST_CHECK_MESSAGE(c2==c_result2,std::string("Inner/Element-Wise Product 1 for ")+typeid(_data_type).name());

}

BOOST_AUTO_TEST_CASE( DenseMIAMultTests )
{

    mult_work<double>(4,3);
    mult_work<float>(4,3);
    mult_work<int>(4,3);
    mult_work<long>(4,3);


}