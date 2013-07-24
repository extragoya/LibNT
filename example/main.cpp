/*Example file on how to get started using LibMIA*/

#include <iostream>



#include "DenseMIA.h"
#include "Index.h"



int main(){


	//declare a set of indices
    LibMIA::MIAINDEX i;
    LibMIA::MIAINDEX j;
    LibMIA::MIAINDEX k;
    LibMIA::MIAINDEX l;
    LibMIA::MIAINDEX m;
    LibMIA::MIAINDEX n;

    //declare our operands
	LibMIA::DenseMIA<double,4> a(5,4,5,4);
	LibMIA::DenseMIA<double,4> b(4,4,5,5);
    //and our resulting MIA
	LibMIA::DenseMIA<double,5> c;


    //initialize DenseMIAs with ones
	a.ones();
    b.ones();

	//calculate a mixed product of inner, outer, and element-wise products
    c(i,k,m,n,l)=a(i,j,k,!l)*b(j,!l,m,n);

}


