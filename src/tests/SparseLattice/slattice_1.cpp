#define EIGEN_SUPERLU_SUPPORT
#include <string>
#include <iostream>
#include "SparseLattice.h"
#include "MappedSparseLattice.h"
#include "Lattice.h"
#include "DenseLattice.h"
#include "slattice_test.h"
typedef LibMIA::SparseLattice<double> sLattice;
typedef LibMIA::MappedSparseLattice<double> mLattice;

int main(){



    dLattice test1= dLattice(4,4,4);
    test1.eye();
    dLattice test2= dLattice(4,4,4);
    test2.ones();

    //dLattice test3=test1.solve(test2);
    //test3.print();

    sLattice stest1(test1);
    //stest1.print();

    sLattice stest2(test2);
    //stest2.print();
    sLattice stest3=stest1.solve(stest2);
    //stest3.print();
//
//    //sLattice stest4=stest1*stest2;
    //
//    stest2.clone(test2);
//    sLattice stest3=stest1*stest2;
//    stest1+=stest2;
//    std::cout <<"********Dense********\n";
//    test3.print();
//    std::cout <<"********Sparse********\n";
//    std::cout <<"Sparse depth" << stest3.depth() <<"\n";
//    stest3.print();



}
