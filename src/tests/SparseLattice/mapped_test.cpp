


#include "slattice_test.h"
////
typedef LibMIA::SparseLattice<double> sLattice;
typedef LibMIA::MappedSparseLattice<double> mLattice;
void mapped_test(){


    dLattice test1= dLattice(4,4,4);
    test1.eye();
    dLattice test2= dLattice(4,4,4);
    test2.ones();

    //dLattice test3=test1.solve(test2);
    //test3.print();

    sLattice stest1(test1);
    stest1.print();

    sLattice stest2(test2);
    stest2.print();
    double * data_1=new double[stest1.size()];
    double * copy_it=data_1;
    long long * index_1=new long long [stest1.size()];
    long long * copy_it_index=index_1;
    std::copy(stest1.data_begin(),stest1.data_end(),&data_1[0]);
    std::copy(stest1.index_begin(),stest1.index_end(),&index_1[0]);
    mLattice mtest1(&data_1[0],&index_1[0],stest1.size(),stest1.height(),stest1.width(),stest1.depth());
    //mtest1.print();

//    double * data_3=new double[stest1.size()];
//    long long * index_3=new long long [stest1.size()];
//    std::copy(stest1.data_begin(),stest1.data_end(),data_3);
//    std::copy(stest1.index_begin(),stest1.index_end(),index_3);
//    mLattice mtest3(data_3,index_3,stest1.size(),stest1.height(),stest1.width(),stest1.depth());
//    mtest3.print();

    double * data_2=new double[stest2.size()];
    long long * index_2=new long long [stest2.size()];
    std::copy(stest2.data_begin(),stest2.data_end(),data_2);
    std::copy(stest2.index_begin(),stest2.index_end(),index_2);
    mLattice mtest2(data_2,index_2,stest2.size(),stest2.height(),stest2.width(),stest2.depth());
    //mtest2.print();
//
    std::cout<<"Got here\n";
    sLattice stest3=mtest1.solve(mtest2);
    stest3.print();

}

