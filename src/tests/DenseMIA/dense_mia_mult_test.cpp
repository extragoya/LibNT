#include <iostream>

#define BOOST_TEST_MODULE DenseMIAMultTests
#include <boost/test/included/unit_test.hpp>

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

    LibMIA::DenseMIA<double,4> a(dim1,dim2,dim1,dim2);
    LibMIA::DenseMIA<double,4> b(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<double,4> c;
    LibMIA::DenseMIA<double,4> c_result(dim1,dim1,dim1,dim1);

    b.ones();

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
    BOOST_CHECK_MESSAGE(c==c_result,"Inner/Outer Product 1");


}

BOOST_AUTO_TEST_CASE( DenseMIAMultTests )
{

    mult_work<double>(4,3);
    mult_work<float>(4,3);
    mult_work<int>(4,3);
    mult_work<long>(4,3);


}
