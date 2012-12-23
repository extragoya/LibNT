#include <string>
#include <iostream>
#include <Lattice.h>
//#include "SparseLattice.h"
#include "DenseLattice.h"
#include "LatticeException.h"

//typedef LibMIA::Lattice<dLattice> Lattice;



int main()
{

    int i=4;

    std::cin >> i;
    //try{
    dLattice test1= dLattice(i,i,i);
    test1.randu(-5,5);
    dLattice test2= dLattice(5,5,5);
    test2.randu(-5,5);
    try{
        dLattice test3=test1.solve(test2);
        test3.print();
    }
    catch(LibMIA::LatticeParameterException& e){
        std::cout << e.what() <<std::endl;

    }

    return 0;
}
