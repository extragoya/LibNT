#include <iostream>

#include "DenseMIA.h"
#include "Index.h"
#include "MIA_Expr.h"

void counter_test();
typedef LibMIA::DenseMIA<double,3> dmia;
typedef LibMIA::DenseMIA<double,4> dmia2;
int main(){

    LibMIA::PRODINDEX i;
    LibMIA::PRODINDEX j;
    LibMIA::PRODINDEX k;
    LibMIA::PRODINDEX l;
    int test=5;

    dmia a(5,7,4);

    dmia2 b(test,7,4,5);
    std::cout<<"HAroo"<<a.at(2,2,2)<<"\n";

    a(i,!j,k)*b(i,!j,k,l);
    counter_test();


}
