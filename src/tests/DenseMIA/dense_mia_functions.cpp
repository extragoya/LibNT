#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE DenseMIAFunctionTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "DenseMIA.h"
#include "Index.h"

template<class _data_type>
void functions_work(){

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;

//    size_t dim1=10;
//    size_t dim2=11;
//    size_t dim3=20;
//    size_t dim4=15;
    size_t dim1=2;
    size_t dim2=3;
    size_t dim3=4;
    size_t dim4=15;

    LibMIA::DenseMIA<_data_type,2> e(dim1,dim2);
    LibMIA::DenseMIA<_data_type,2> e2(dim1,dim2);
    LibMIA::DenseMIA<_data_type,2> f;
    LibMIA::DenseMIA<_data_type,3> a(dim1,dim2,dim3);
    LibMIA::DenseMIA<_data_type,3> a2(dim1,dim2,dim3);
    LibMIA::DenseMIA<_data_type,3> b;
    LibMIA::DenseMIA<_data_type,4> c(dim1,dim2,dim3,dim4);
    LibMIA::DenseMIA<_data_type,4> c2(dim1,dim2,dim3,dim4);
    LibMIA::DenseMIA<_data_type,4> d;

    //****************In place permutation
    e.randu(-50,50);
    e2=e;
    f(j,i)=e(i,j);
    e2.inplace_permute(1,0);
    //e.print();
    //f.print();
    //e2.print();
    BOOST_CHECK_MESSAGE(f==e2,std::string("Second-order in place permutation 1 for ")+typeid(_data_type).name());

    a.randu(-50,50);
    //a.print();
    a2=a;

    b(j,k,i)=a(i,j,k);

    a2.inplace_permute(1,2,0);

    //b.print();
    //a2.print();
    BOOST_CHECK_MESSAGE(b==a2,std::string("Third-order in place permutation 1 for ")+typeid(_data_type).name());
    a2=a;
    b(k,i,j)=a(i,j,k);
    a2.inplace_permute(2,0,1);
    BOOST_CHECK_MESSAGE(b==a2,std::string("Third-order in place permutation 2 for ")+typeid(_data_type).name());

    c.randu(-50,50);
    c2=c;
    d(k,l,i,j)=c(i,j,k,l);
    c2.inplace_permute(2,3,0,1);
    BOOST_CHECK_MESSAGE(d==c2,std::string("Fourth-order in place permutation 1 for ")+typeid(_data_type).name());

    c2=c;
    d(j,l,i,k)=c(i,j,k,l);
    c2.inplace_permute(1,3,0,2);
    BOOST_CHECK_MESSAGE(d==c2,std::string("Fourth-order in place permutation 2 for ")+typeid(_data_type).name());

    auto lat1=d.toLatticeCopy(std::array<size_t, 1>{{3}},std::array<size_t, 2>{{2,1}},std::array<size_t, 1>{{0}});
    auto lat2=d.toLatticePermute(std::array<size_t, 1>{{3}},std::array<size_t, 2>{{2,1}},std::array<size_t, 1>{{0}});
    BOOST_CHECK_MESSAGE(lat1==lat2,std::string("Lattice Mapping test for ")+typeid(_data_type).name());

}

BOOST_AUTO_TEST_CASE( DenseMIAFunctionTests )
{



    functions_work<double>();

//    functions_work<float>();
//    functions_work<int>();
//    functions_work<long>();






}
