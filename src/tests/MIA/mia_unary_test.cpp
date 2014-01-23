#include <iostream>

#define BOOST_TEST_MODULE MIAUnaryTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"
#include "SparseMIA.h"
#include "LibMIAHelpers.h"
#include "Index.h"

//typedef LibMIA::DenseMIA<double,3> dmia;

template<class data_type>
void sparse_unary_work(){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;
    size_t _dim=5;

    LibMIA::SparseMIA<data_type,6> a(_dim,_dim+1,_dim,_dim+1,_dim-1,_dim-2);

    LibMIA::SparseMIA<data_type,4> b;
    LibMIA::SparseMIA<data_type,4> other_b;
    LibMIA::SparseMIA<data_type,3> c;
    LibMIA::SparseMIA<data_type,3> other_c;
    LibMIA::SparseMIA<data_type,2> d;
    LibMIA::SparseMIA<data_type,2> other_d;
    LibMIA::SparseMIA<data_type,1> e;
    LibMIA::SparseMIA<data_type,1> other_e;
    LibMIA::SparseMIA<data_type,5> f;
    LibMIA::SparseMIA<data_type,5> other_f;
    auto delta=LibMIA::create_delta<data_type,2>(_dim);
    auto delta2=LibMIA::create_delta<data_type,2>(_dim+1);
    auto delta_3=LibMIA::create_delta<data_type,3>(_dim);

    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();



    b(j,k,l,m)=a(i,j,i,k,l,m); //contraction
    other_b(j,k,l,m)=a(i,j,n,k,l,m)*delta(i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Simple Contraction 1 for Sparse ")+typeid(data_type).name());

    c=LibMIA::SparseMIA<data_type,3>(_dim,_dim+1,_dim);
    c.resize(c.dimensionality()/3);
    c.randu(-5,5);
    c.rand_indices();
    c.collect_duplicates();
    e(i)=c(j,i,j);
    other_e(i)=c(j,i,k)*delta(j,k);
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Simple Contraction 2 for Sparse ")+typeid(data_type).name());


//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+2,_dim-2);
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    c(j,l,m)=a(i,j,i,i,l,m); //contraction
    other_c(j,l,m)=a(i,j,n,k,l,m)*delta_3(i,n,k); //equivalent operation
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(other_c,test_precision<data_type>()),std::string("Triple Contraction 1 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    d(l,m)=a(i,j,i,j,l,m); //contraction
    other_d(l,m)=a(i,j,n,k,l,m)*delta(i,n)*delta2(j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Complex Contraction 1 for Sparse ")+typeid(data_type).name());

    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    e(m)=a(i,j,i,j,i,m); //contraction
    other_e(m)=a(i,j,n,k,l,m)*delta_3(i,n,l)*delta2(j,k); //equivalent operation
//    e.print();
//    other_e.print();
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Complex Contraction 2 for Sparse ")+typeid(data_type).name());

    d=LibMIA::SparseMIA<data_type,2>(_dim,_dim);
    d.resize(a.dimensionality()/3);
    d.randu(-5,5);
    d.rand_indices();
    d.collect_duplicates();
    e(i)=d(!i,!i);
    other_e(i)=d(!i,j)*delta(!i,j);
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Very Simple Attraction 1 for Sparse ")+typeid(data_type).name());
//
    c=LibMIA::SparseMIA<data_type,3>(_dim,_dim+1,_dim); //now try with non-uniform dimensions
    c.resize(a.dimensionality()/3);
    c.randu(-5,5);
    c.rand_indices();
    c.collect_duplicates();
    d(i,j)=c(!i,j,!i); //attraction
    other_d(i,j)=c(!i,j,n)*delta(!i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Simple Attraction 1 for Sparse ")+typeid(data_type).name());
//
//
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+2,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    f(i,j,k,l,m)=a(!i,j,!i,k,l,m); //attraction
    other_f(i,j,k,l,m)=a(!i,j,n,k,l,m)*delta(!i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(f.fuzzy_equals(other_f,test_precision<data_type>()),std::string("Simple Attraction 2 for Sparse ")+typeid(data_type).name());

//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+2,_dim-2);
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    b(i,j,l,m)=a(!i,j,!i,!i,l,m); //attraction
    other_b(i,j,l,m)=a(!i,j,n,k,l,m)*delta_3(!i,n,k); //equivalent operation
    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Triple Attraction 1 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    b(i,j,l,m)=a(!i,!j,!i,!j,l,m); //attraction
    other_b(i,j,l,m)=a(!i,!j,n,k,l,m)*delta(!i,n)*delta2(!j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Complex Attraction 1 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    c(i,j,m)=a(!i,!j,!i,!j,!i,m); //attraction
    other_c(i,j,m)=a(!i,!j,n,k,l,m)*delta_3(!i,n,l)*delta2(!j,k); //equivalent operation
    //c.print();
    //other_c.print();
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(other_c,test_precision<data_type>()),std::string("Complex Attraction 2 for Sparse ")+typeid(data_type).name());
//
    b=LibMIA::SparseMIA<data_type,4>(_dim,_dim+1,_dim,_dim+1); //now try with non-uniform dimensions
    b.resize(a.dimensionality()/3);
    b.randu(-5,5);
    b.rand_indices();
    b.collect_duplicates();
    d(i,j)=b(!i,!j,!i,!j); //attraction
    other_d(i,j)=b(!i,!j,n,k)*delta(!i,n)*delta2(!j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Complex Attraction 3 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    c(i,k,m)=a(!i,j,!i,j,k,m); //attraction
    other_c(i,k,m)=a(n,l,!i,j,k,m)*delta(!i,n)*delta2(j,l); //attraction
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(other_c,test_precision<data_type>()),std::string("Combined Contraction/Attraction 1 for Sparse ")+typeid(data_type).name());


    b=LibMIA::SparseMIA<data_type,4>(_dim,_dim+1,_dim,_dim+1); //now try with non-uniform dimensions
    b.resize(a.dimensionality()/3);
    b.randu(-5,5);
    b.rand_indices();
    b.collect_duplicates();
    e(i)=b(!i,j,!i,j); //attraction
    other_e(i)=b(n,l,!i,j)*delta(!i,n)*delta2(j,l); //attraction
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Combined Contraction/Attraction 2 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+1,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    d(i,m)=a(!i,j,!i,!i,j,m); //attraction
    other_d(i,m)=a(!i,j,n,k,l,m)*delta_3(!i,n,k)*delta2(j,l); //attraction
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Combined Contraction/Attraction 3 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+1,_dim-2); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    d(j,m)=a(i,!j,i,i,!j,m); //attraction
    other_d(j,m)=a(i,!j,n,k,l,m)*delta_3(i,n,k)*delta2(!j,l); //attraction
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Combined Contraction/Attraction 4 for Sparse ")+typeid(data_type).name());
//
    a=LibMIA::SparseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim); //now try with non-uniform dimensions
    a.resize(a.dimensionality()/3);
    a.randu(-5,5);
    a.rand_indices();
    a.collect_duplicates();
    d(i,k)=a(!i,j,!i,j,!k,!k); //attraction
    other_d(i,k)=a(n,l,!i,j,!k,m)*delta(!i,n)*delta2(j,l)*delta(!k,m); //attraction
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Combined Contraction/Attraction 5 for Sparse ")+typeid(data_type).name());



//

}


template<class data_type>
void dense_unary_work(){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;
    size_t _dim=5;

    LibMIA::DenseMIA<data_type,6> a(_dim,_dim,_dim,_dim,_dim,_dim);

    LibMIA::DenseMIA<data_type,4> b;
    LibMIA::DenseMIA<data_type,4> other_b;
    LibMIA::DenseMIA<data_type,3> c;
    LibMIA::DenseMIA<data_type,3> other_c;
    LibMIA::DenseMIA<data_type,2> d;
    LibMIA::DenseMIA<data_type,2> other_d;
    LibMIA::DenseMIA<data_type,1> e;
    LibMIA::DenseMIA<data_type,1> other_e;
    LibMIA::DenseMIA<data_type,5> f;
    LibMIA::DenseMIA<data_type,5> other_f;
    auto delta=LibMIA::create_delta<data_type,2>(_dim);
    auto delta2=LibMIA::create_delta<data_type,2>(_dim+1);
    auto delta_3=LibMIA::create_delta<data_type,3>(_dim);
    a.randu(-5,5);


    const LibMIA::DenseMIA<data_type,6> temp_a(a); //check const_correctness
    b(j,k,l,m)=temp_a(i,j,i,k,l,m); //contraction
    other_b(j,k,l,m)=temp_a(i,j,n,k,l,m)*delta(i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Simple Contraction 1 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+2,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    b(j,k,l,m)=a(i,j,i,k,l,m); //contraction
    other_b(j,k,l,m)=a(i,j,n,k,l,m)*delta(i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Simple Contraction 2 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+2,_dim-2);
    a.randu(-5,5);
    c(j,l,m)=a(i,j,i,i,l,m); //contraction
    other_c(j,l,m)=a(i,j,n,k,l,m)*delta_3(i,n,k); //equivalent operation
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(other_c,test_precision<data_type>()),std::string("Triple Contraction 1 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    d(l,m)=a(i,j,i,j,l,m); //contraction
    other_d(l,m)=a(i,j,n,k,l,m)*delta(i,n)*delta2(j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Complex Contraction 1 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    e(m)=a(i,j,i,j,i,m); //contraction
    other_e(m)=a(i,j,n,k,l,m)*delta_3(i,n,l)*delta2(j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Complex Contraction 2 for Dense ")+typeid(data_type).name());

    d=LibMIA::DenseMIA<data_type,2>(_dim,_dim);
    d.randu(-5,5);
    e(i)=d(!i,!i);
    other_e(i)=d(!i,j)*delta(!i,j);
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Very Simple Attraction 1 for Dense ")+typeid(data_type).name());

    c=LibMIA::DenseMIA<data_type,3>(_dim,_dim+1,_dim); //now try with non-uniform dimensions
    c.randu(-5,5);
    d(i,j)=c(!i,j,!i); //attraction
    other_d(i,j)=c(!i,j,n)*delta(!i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Simple Attraction 1 for Dense ")+typeid(data_type).name());



    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+2,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    f(i,j,k,l,m)=a(!i,j,!i,k,l,m); //attraction
    other_f(i,j,k,l,m)=a(!i,j,n,k,l,m)*delta(!i,n); //equivalent operation
    BOOST_CHECK_MESSAGE(f.fuzzy_equals(other_f,test_precision<data_type>()),std::string("Simple Attraction 2 for Dense ")+typeid(data_type).name());


    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+2,_dim-2);
    a.randu(-5,5);
    b(i,j,l,m)=a(!i,j,!i,!i,l,m); //attraction
    other_b(i,j,l,m)=a(!i,j,n,k,l,m)*delta_3(!i,n,k); //equivalent operation
    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Triple Attraction 1 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim-1,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    b(i,j,l,m)=a(!i,!j,!i,!j,l,m); //attraction
    other_b(i,j,l,m)=a(!i,!j,n,k,l,m)*delta(!i,n)*delta2(!j,k); //equivalent operation

    BOOST_CHECK_MESSAGE(b.fuzzy_equals(other_b,test_precision<data_type>()),std::string("Complex Attraction 1 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    c(i,j,m)=a(!i,!j,!i,!j,!i,m); //attraction
    other_c(i,j,m)=a(!i,!j,n,k,l,m)*delta_3(!i,n,l)*delta2(!j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(other_c,test_precision<data_type>()),std::string("Complex Attraction 2 for Dense ")+typeid(data_type).name());

    b=LibMIA::DenseMIA<data_type,4>(_dim,_dim+1,_dim,_dim+1); //now try with non-uniform dimensions
    b.randu(-5,5);
    d(i,j)=b(!i,!j,!i,!j); //attraction
    other_d(i,j)=b(!i,!j,n,k)*delta(!i,n)*delta2(!j,k); //equivalent operation
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Complex Attraction 3 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    c(i,k,m)=a(!i,j,!i,j,k,m); //attraction
    other_c(i,k,m)=a(n,l,!i,j,k,m)*delta(!i,n)*delta2(j,l); //attraction
    BOOST_CHECK_MESSAGE(c.fuzzy_equals(other_c,test_precision<data_type>()),std::string("Combined Contraction/Attraction 1 for Dense ")+typeid(data_type).name());

    b=LibMIA::DenseMIA<data_type,4>(_dim,_dim+1,_dim,_dim+1); //now try with non-uniform dimensions
    b.randu(-5,5);
    e(i)=b(!i,j,!i,j); //attraction
    other_e(i)=b(n,l,!i,j)*delta(!i,n)*delta2(j,l); //attraction
    BOOST_CHECK_MESSAGE(e.fuzzy_equals(other_e,test_precision<data_type>()),std::string("Combined Contraction/Attraction 2 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+1,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    d(i,m)=a(!i,j,!i,!i,j,m); //attraction
    other_d(i,m)=a(!i,j,n,k,l,m)*delta_3(!i,n,k)*delta2(j,l); //attraction
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Combined Contraction/Attraction 3 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim,_dim+1,_dim-2); //now try with non-uniform dimensions
    a.randu(-5,5);
    d(j,m)=a(i,!j,i,i,!j,m); //attraction
    other_d(j,m)=a(i,!j,n,k,l,m)*delta_3(i,n,k)*delta2(!j,l); //attraction
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Combined Contraction/Attraction 4 for Dense ")+typeid(data_type).name());

    a=LibMIA::DenseMIA<data_type,6>(_dim,_dim+1,_dim,_dim+1,_dim,_dim); //now try with non-uniform dimensions
    a.randu(-5,5);
    d(i,k)=a(!i,j,!i,j,!k,!k); //attraction
    other_d(i,k)=a(n,l,!i,j,!k,m)*delta(!i,n)*delta2(j,l)*delta(!k,m); //attraction
    BOOST_CHECK_MESSAGE(d.fuzzy_equals(other_d,test_precision<data_type>()),std::string("Combined Contraction/Attraction 5 for Dense ")+typeid(data_type).name());
//    d.print();
//    other_d.print();



//

}

BOOST_AUTO_TEST_CASE( MIAUnaryTests )
{

    dense_unary_work<double>();
    dense_unary_work<float>();
    dense_unary_work<int>();
    dense_unary_work<long>();

    sparse_unary_work<double>();
    sparse_unary_work<float>();
    sparse_unary_work<int>();
    sparse_unary_work<long>();



}
