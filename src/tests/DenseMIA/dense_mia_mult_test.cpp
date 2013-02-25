#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE DenseMIAMultTests
#include "MIAConfig.h"
#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

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

    LibMIA::DenseMIA<_data_type,4> a(dim1,dim2,dim1,dim2);
    LibMIA::DenseMIA<_data_type,4> b(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,2> b2(dim2,dim2);

    LibMIA::DenseMIA<_data_type,4> c;
    LibMIA::DenseMIA<_data_type,2> c2;
    LibMIA::DenseMIA<_data_type,4> c_result(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,2> c_result2(dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> c_result3(dim1,dim2,dim1,dim2);

    LibMIA::DenseMIA<_data_type,1> d(dim2);
    LibMIA::DenseMIA<_data_type,2> d2(dim2,dim2);

    b.ones();

    //set each the values of a to increment along inner product indices
    size_t val,entries=((dim2*dim2+1)*dim2*dim2)/2;
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

    //test inner and element-wise product. Have b be assigned a scalar value that increases while traversing i
    val=0;
    for(size_t _i=0;_i<dim1;++_i){
        val++;
        for(size_t _j=0;_j<dim1;++_j){
            c_result2.at(_i,_j)=val*entries;
            for(size_t _k=0;_k<dim2;++_k){
                for(size_t _l=0;_l<dim2;++_l){
                    b.at(_k,_l,_i,_j)=val;
                }
            }
        }
    }

    c2(i,j)=a(!i,k,!j,l)*b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2==c_result2,std::string("Inner/Element-Wise Product 1 for ")+typeid(_data_type).name());
    c2(i,j)=a(!i,l,!j,k)*b(k,l,!i,!j);
    BOOST_CHECK_MESSAGE(c2==c_result2,std::string("Inner/Element-Wise Product 2 for ")+typeid(_data_type).name());
    c2(j,i)=a(!j,l,!i,k)*b(k,l,!j,!i);
    BOOST_CHECK_MESSAGE(c2==c_result2,std::string("Inner/Element-Wise Product 3 for ")+typeid(_data_type).name());

    //test inter outer products
    val=1;
    a.ones();
    for(size_t _i=0;_i<dim2;++_i){

        for(size_t _j=0;_j<dim2;++_j){
            b2.at(_i,_j)=val++;
        }
    }

    c_result3.zeros();
    for(size_t _i=0;_i<dim1;++_i){
        for(size_t _k=0;_k<dim1;++_k){
            val=1;
            for(size_t _j=0;_j<dim2;++_j){
                for(size_t _l=0;_l<dim2;++_l){
                    c_result3.at(_i,_j,_k,_l)=val++;
                }
            }
        }
    }
    c(i,j,k,l)=a(i,!j,k,!l)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c==c_result3,std::string("Outer/Element-Wise Product 1 for ")+typeid(_data_type).name());
    c(i,j,k,l)=a(k,!j,i,!l)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c==c_result3,std::string("Outer/Element-Wise Product 2 for ")+typeid(_data_type).name());

    c(i,j,k,l)=a(k,!l,i,!j)*b2(!j,!l);
    BOOST_CHECK_MESSAGE(c==c_result3,std::string("Outer/Element-Wise Product 3 for ")+typeid(_data_type).name());

    val=1;
    for(size_t _i=0;_i<dim2;++_i)
        d.at(_i)=val++;

    size_t val2;
    c_result3.zeros();
    for(size_t _i=0;_i<dim1;++_i){
        for(size_t _k=0;_k<dim1;++_k){
            val=1;
            for(size_t _j=0;_j<dim2;++_j){
                val2=1;
                for(size_t _l=0;_l<dim2;++_l){
                    c_result3.at(_i,_j,_k,_l)=val++*val2++;
                }
            }
        }
    }


    c(i,j,k,l)=~(a(i,!j,k,!!l)*b2(!j,!!l))*d(!l);
    BOOST_CHECK_MESSAGE(c==c_result3,std::string("Repeated Element-Wise Product 1 for ")+typeid(_data_type).name());

    //repeat inner/outer product test but use a ternary inner product instead
    d2.ones();
    b.ones();
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
    c(i,k,m,n)=~(a(i,!j,k,!l)*b(!j,!l,m,n))*d2(j,l);
    BOOST_CHECK_MESSAGE(c==c_result,std::string("Ternary Inner Product 1 for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(a(i,!j,k,!l)*b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c==c_result,std::string("Ternary Inner Product 2 for ")+typeid(_data_type).name() );
    c(i,k,m,n)=~(a(i,!l,k,!j)*b(!j,!l,m,n))*d2(l,j);
    BOOST_CHECK_MESSAGE(c==c_result,std::string("Ternary Inner Product 2 for ")+typeid(_data_type).name() );

}

BOOST_AUTO_TEST_CASE( DenseMIAMultTests )
{

    mult_work<double>(3,3);
    mult_work<float>(3,3);
    mult_work<int>(3,3);
    mult_work<long>(3,3);


    mult_work<double>(4,3);
    mult_work<float>(4,3);
    mult_work<int>(4,3);
    mult_work<long>(4,3);


    mult_work<double>(3,4);
    mult_work<float>(3,4);
    mult_work<int>(3,4);
    mult_work<long>(3,4);


}
