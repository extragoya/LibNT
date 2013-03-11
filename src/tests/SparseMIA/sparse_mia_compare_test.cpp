#include <string>

#define BOOST_TEST_MODULE SparseCompareTests




#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include "DenseMIA.h"
#include "SparseMIA.h"






template<typename MIAType, typename DenseMIAType>
void do_work(size_t dim1, size_t dim2, size_t dim3){


//    MIA_Type test1(dim1,dim2,dim3);
//    test1.resize(test1.dimensionality()/2);




}




BOOST_AUTO_TEST_CASE( SparseCompareTests )
{

    size_t dim1=3,dim2=4,dim3=5;
    //testing to make sure a copy constructor doesn't get called erroneously
    do_work<LibMIA::SparseMIA<float,3>, LibMIA::DenseMIA<float,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<double,3>, LibMIA::DenseMIA<double,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<int32_t,3>,LibMIA::DenseMIA<int32_t,3> >(dim1,dim2,dim3);
    do_work<LibMIA::SparseMIA<int64_t,3>,LibMIA::DenseMIA<int64_t,3> >(dim1,dim2,dim3);
}
