#include <string>

#define BOOST_TEST_MODULE MixedMIAMergeTests




#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"
#include "SparseMIA.h"
constexpr int dim=5;



template<typename data_type>
void do_work(size_t dim1){

    LibMIA::DenseMIA<data_type,3> dense_a(dim1,dim1,dim1);
    LibMIA::DenseMIA<data_type,3> dense_b(dim1,dim1,dim1);
    LibMIA::DenseMIA<data_type,3> dense_c;
    LibMIA::DenseMIA<data_type,3> dense_c2;

    LibMIA::SparseMIA<data_type,3> sparse_a;
    LibMIA::SparseMIA<data_type,3> sparse_c2;

    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;

    dense_a.randu(0,20);
    dense_b.randu(0,20);
    for(auto it=dense_a.data_begin();it<dense_a.data_end();++it)
        if(*it<15)
            *it=0;
    for(auto it=dense_b.data_begin();it<dense_b.data_end();++it)
        if(*it<15)
            *it=0;

    sparse_a=dense_a;
    dense_c(i,j,k)=dense_a(i,j,k)+dense_b(i,j,k);
    dense_c2(i,j,k)=sparse_a(i,j,k)+dense_b(i,j,k);

    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Add 1a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c2(i,j,k)=dense_b(i,j,k)+sparse_a(i,j,k);
    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Add 1b for ")+typeid(data_type).name());

    dense_c(k,i,j)=dense_a(i,k,j)+dense_b(j,k,i);
    dense_c2(k,i,j)=sparse_a(i,k,j)+dense_b(j,k,i);

    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Add 2a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c2(k,i,j)=dense_b(j,k,i)+sparse_a(i,k,j);
    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Add 2b for ")+typeid(data_type).name());

    //**Destructive
    dense_c=dense_a;
    sparse_c2=sparse_a;
    dense_c(i,j,k)+=dense_b(i,j,k);
    sparse_c2(i,j,k)+=dense_b(i,j,k);

    BOOST_CHECK_MESSAGE(dense_c==sparse_c2,std::string("Destructive Add 1a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c=dense_b;
    dense_c2=dense_b;
    dense_c(i,j,k)+=dense_a(i,j,k);
    dense_c2(i,j,k)+=sparse_a(i,j,k);
    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Destructive Add 1b for ")+typeid(data_type).name());


    dense_c=dense_a;
    sparse_c2=sparse_a;
    dense_c(k,i,j)+=dense_b(j,k,i);
    sparse_c2(k,i,j)+=dense_b(j,k,i);

    BOOST_CHECK_MESSAGE(dense_c==sparse_c2,std::string("Destructive Add 2a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c=dense_b;
    dense_c2=dense_b;
    dense_c(k,i,j)+=dense_a(i,k,j);
    dense_c2(k,i,j)+=sparse_a(i,k,j);

    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Destructive Add 2b for ")+typeid(data_type).name());

    //******Subtract******************

    dense_c(i,j,k)=dense_a(i,j,k)-dense_b(i,j,k);
    dense_c2(i,j,k)=sparse_a(i,j,k)-dense_b(i,j,k);

    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Sub 1a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c(i,j,k)=dense_b(i,j,k)-dense_a(i,j,k);
    dense_c2(i,j,k)=dense_b(i,j,k)-sparse_a(i,j,k);
    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Sub 1b for ")+typeid(data_type).name());

    dense_c(k,i,j)=dense_a(i,k,j)-dense_b(j,k,i);
    dense_c2(k,i,j)=sparse_a(i,k,j)-dense_b(j,k,i);

    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Sub 2a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c(k,i,j)=dense_b(j,k,i)-dense_a(i,k,j);
    dense_c2(k,i,j)=dense_b(j,k,i)-sparse_a(i,k,j);
    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Non-destructive Sub 2b for ")+typeid(data_type).name());


    //**Destructive
    dense_c=dense_a;
    sparse_c2=sparse_a;
    dense_c(i,j,k)-=dense_b(i,j,k);
    sparse_c2(i,j,k)-=dense_b(i,j,k);

    BOOST_CHECK_MESSAGE(dense_c==sparse_c2,std::string("Destructive Sub 1a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c=dense_b;
    dense_c2=dense_b;
    dense_c(i,j,k)-=dense_a(i,j,k);
    dense_c2(i,j,k)-=sparse_a(i,j,k);
    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Destructive Sub 1b for ")+typeid(data_type).name());


    dense_c=dense_a;
    sparse_c2=sparse_a;
    dense_c(k,i,j)-=dense_b(j,k,i);
    sparse_c2(k,i,j)-=dense_b(j,k,i);

    BOOST_CHECK_MESSAGE(dense_c==sparse_c2,std::string("Destructive Sub 2a for ")+typeid(data_type).name());
    //switch order of operands
    dense_c=dense_b;
    dense_c2=dense_b;
    dense_c(k,i,j)-=dense_a(i,k,j);
    dense_c2(k,i,j)-=sparse_a(i,k,j);

    BOOST_CHECK_MESSAGE(dense_c==dense_c2,std::string("Destructive Sub 2b for ")+typeid(data_type).name());


}





BOOST_AUTO_TEST_CASE( MixedMIAMergeTests )
{

    size_t dim1=3;


    do_work<float>(dim1);
    do_work<double>(dim1);
    do_work<int32_t>(dim1);
    do_work<int64_t>(dim1);

}



