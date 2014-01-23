#include <iostream>
#include <typeinfo>

#include "MIAConfig.h"
#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "SparseMIA.h"
#include "DenseMIA.h"
#include "Index.h"
#include "LibMIAException.h"
#include "LibMIAUtil.h"

template<typename data_type, size_t _order>
void random_mia(LibMIA::DenseMIA<data_type,_order> & mia,double _prob,bool need_ranked=true)
{
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(LibMIA::LibMIA_gen(), uni_dist);
    boost::uniform_real<> uni_dist2(-10,10);
    boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni2(LibMIA::LibMIA_gen(), uni_dist2);
    mia.zeros();
    for(auto it =mia.data_begin(); it<mia.data_end(); ++it)
    {

        if(uni()<_prob)
        {
            *it=uni2();
        }

    }
}


template<class _data_type>
void solve_work(size_t dim1, size_t dim2){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;
    //LibMIA::MIAINDEX o;
    //LibMIA::MIAINDEX p;

    LibMIA::DenseMIA<_data_type,4> dense_a(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_a2(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_b(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_b2(dim2,dim2,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_c(dim1,dim1,dim1,dim1);
    LibMIA::DenseMIA<_data_type,4> dense_d2(dim1,dim1,dim1,dim1);

    LibMIA::SparseMIA<_data_type,4> a(dim1,dim1,dim1,dim1);
    LibMIA::SparseMIA<_data_type,4> a2(dim2,dim2,dim1,dim1);
    LibMIA::SparseMIA<_data_type,4> b(dim1,dim1,dim1,dim1);
    LibMIA::SparseMIA<_data_type,4> b2(dim2,dim2,dim1,dim1);
    //result of a solution of sparse equations is set to dense
    LibMIA::DenseMIA<_data_type,4> c(dim1,dim1,dim1,dim1);


    double _prob=0.4;
    bool flag=true;
    random_mia(dense_b,_prob);
    while(flag){
        random_mia(dense_a,_prob);

        dense_c(i,j,m,n)=dense_a(i,j,k,l)|dense_b(k,l,m,n);
        if (dense_c.solveInfo()==LibMIA::FullyRanked)
            flag=false;

    }
    const LibMIA::DenseMIA<_data_type,4> temp_a(dense_a); //make sure we've sorted out const correctness
    const LibMIA::DenseMIA<_data_type,4> temp_b(dense_b);
    a=dense_a;
    b=dense_b;
    c(i,j,m,n)=temp_a(i,j,k,l)|b(k,l,m,n);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product Inverse 1 (dense) for ")+typeid(_data_type).name() );

    c(i,j,m,n)=a(i,j,k,l)|temp_b(k,l,m,n);
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product Inverse 1 (sparse) for ")+typeid(_data_type).name() );

    flag=true;
    while(flag){


        dense_c(i,j,m,n)=dense_a(i,k,j,l)|dense_b(k,l,m,n);
        if (dense_c.solveInfo()==LibMIA::FullyRanked)
            break;
        random_mia(dense_a,_prob);
    }
    a=dense_a;
    c(i,j,m,n)=dense_a(i,k,j,l)|b(k,l,m,n);

    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product Inverse 2 (dense) for ")+typeid(_data_type).name() );

    c(i,j,m,n)=a(i,k,j,l)|dense_b(k,l,m,n);

    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer Product Inverse 2 (sparse) for ")+typeid(_data_type).name() );

    _prob=0.7;
    random_mia(dense_a,_prob);
    flag=true;
    while(flag){


        dense_c(i,j,l,m)=dense_a(i,!j,k,!l)|dense_b(k,!j,!l,m);
        if (dense_c.solveInfo()==LibMIA::FullyRanked)
            break;

        random_mia(dense_a,_prob);
    }
    a=dense_a;
    c(i,j,l,m)=dense_a(i,!j,k,!l)|b(k,!j,!l,m);

    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer/Inter Product Inverse 1 (dense) for ")+typeid(_data_type).name() );

    c(i,j,l,m)=a(i,!j,k,!l)|dense_b(k,!j,!l,m);

    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer/Inter Product Inverse 1 (sparse) for ")+typeid(_data_type).name() );

    flag=true;
    while(flag){


        dense_c(i,j,l,m)=dense_a(!i,!j,k,l)|dense_b(m,!j,k,!i);
        if (dense_c.solveInfo()==LibMIA::FullyRanked)
            break;
        random_mia(dense_a,_prob);
    }
    a=dense_a;
    c(i,j,l,m)=dense_a(!i,!j,k,l)|b(m,!j,k,!i);

    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer/Inter Product Inverse 2 (dense) for ")+typeid(_data_type).name() );

    c(i,j,l,m)=a(!i,!j,k,l)|dense_b(m,!j,k,!i);

    BOOST_CHECK_MESSAGE(c.fuzzy_equals(dense_c,test_precision<_data_type>()),std::string("Inner/Outer/Inter Product Inverse 2 (sparse) for ")+typeid(_data_type).name() );



//    c(i,j,m,n)=a2(k,l,i,j)|b2(k,l,m,n);
//    //test with normal equations
//    d(o,p,m,n)=a2(k,l,o,p)*a2(k,l,i,j)*c(i,j,m,n);
//    d2(i,j,m,n)=a2(k,l,i,j)*b2(k,l,m,n);
//    BOOST_CHECK_MESSAGE(d.fuzzy_equals(d2,test_precision<_data_type>()),std::string("Inner/Outer Product Least Squares 1 for ")+typeid(_data_type).name() );

}


