#define BOOST_TEST_MODULE SparseLatticeSolveTestsFLOAT








#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "sparse_lattice_solve_test.hpp"


BOOST_AUTO_TEST_CASE( SparseLatticeSolveTestsFLOAT )
{


    //multwork<double>(3,3,3,3);

    //solvework<double>(10,10,10,1);

    //solvework<double>(20,20,10,20);

    solvework<float>(20,20,10,20);
    //solvework<float>(5,5,4,5);




}
