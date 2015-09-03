/*Example file on how to get started using LibNT*/

#include <iostream>


#include "SparseMIA.h"
#include "DenseMIA.h"
#include "Index.h"



int main(){

    using namespace LibMIA;
	//declare a set of indices (must use the MIAINDEX macro, line by line)
    NTINDEX i;
	NTINDEX j;
	NTINDEX k;
	NTINDEX l;
	NTINDEX m;
	NTINDEX n;

    //declare our operands along with their size. Template arguments are datatype, then MIA degree (or order)
	DenseNT<double,4> a(5,5,5,5);
	DenseNT<double, 4> b(5, 5, 5, 5);
    //and our resulting MIA
	DenseNT<double, 5> c;
	DenseNT<double, 4> c2;

    //initialize DenseMIAs with random values
	a.randu(-2,2);
    b.randu(-2,2);

    //modify a value within an array
    a.at(1,2,1,2)=2;

	//calculate a mixed product of inner, outer, and element-wise products
    c(i,k,m,n,l)=a(i,j,k,!l)*b(j,!l,m,n);

    //perform a solution of equations
    c2(i,j,k,l)=a(m,n,i,j)|b(m,n,k,l);

    SparseNT<double,3> d(4,4,4), e(4,4,4);

    //push back data, index pairs (need to do some work on making this easier, as right now it's a linearized index)
    d.push_back(5,0);
    d.push_back(3,20);
    d.push_back(-6,5);

    e.push_back(8,8);
    e.push_back(-4,30);
    e.push_back(10.3,17);

    //perform a destructive add
    e(i,j,k)+=d(k,i,j);
    e.print();

}


