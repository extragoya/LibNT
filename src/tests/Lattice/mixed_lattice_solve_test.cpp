#define BOOST_TEST_MODULE MixedLatticeSolveTests



#include "DenseLattice.h"
#include "SparseLattice.h"
#include "MIAConfig.h"
#include "LibMIAUtil.h"

#include <Eigen/Dense>

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif


template<typename data_type>
void random_matrix(LibMIA::DenseLattice<data_type> & lat,double _prob,bool need_ranked=true, bool LSQR=false)
{
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(LibMIA::gen, uni_dist);
    boost::uniform_real<> uni_dist2(-10,10);
    boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni2(LibMIA::gen, uni_dist2);
    lat.zeros();
    for(size_t k=0;k<(size_t)lat.depth();++k){
        auto _start=lat.data_begin()+k*lat.width()*lat.height();
        auto _end=lat.data_begin()+(k+1)*lat.width()*lat.height();
        bool flag=true;
        while(flag){
            for(auto it=_start;it<_end;++it){
                if(uni()<_prob){
                    *it=uni2();
                }

            }
            if(need_ranked){
                if(LSQR){
                    auto _SVD=lat.tab_matrix(k).jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
                    if(!_SVD.nonzeroSingularValues()!=lat.width())
                        flag=false;
                }
                else{
                    auto _QR=lat.tab_matrix(k).colPivHouseholderQr();
                    if(_QR.rank()==lat.width())
                        flag=false;
                }
            }
            else
                flag=false;



        }
    }
}

template<typename data_type>
void solvework(size_t m1, size_t n1, size_t n2, size_t p){

    typedef LibMIA::DenseLattice<data_type> denseType;
    typedef LibMIA::SparseLattice<data_type> sparseType;


    denseType DenseLat1(m1,n1,p);
    denseType DenseLat1_lsqr(m1,2*n1,p);
    denseType DenseLat2(n1,n2,p);

    denseType DenseLat3;
    denseType DenseLat3_mixed;

    sparseType SparseLat1(m1,n1,p);
    sparseType SparseLat2(n1,n2,p);
    sparseType SparseLat3;

    random_matrix(DenseLat1,0.4);
    random_matrix(DenseLat1_lsqr,0.4,false,true);
    random_matrix(DenseLat2,0.4,false);

    SparseLat1=DenseLat1;
    SparseLat2=DenseLat2;

    DenseLat3=DenseLat1.solve(DenseLat2);
    DenseLat3_mixed=DenseLat1.solve(SparseLat2);



    BOOST_CHECK_MESSAGE(DenseLat3.fuzzy_equals(DenseLat3_mixed,test_precision<data_type>()),std::string("Full Dimension Solve Test 1 for ")+typeid(data_type).name());
    //now invert the sparse version instead
    DenseLat3_mixed=SparseLat1.solve(DenseLat2);
    BOOST_CHECK_MESSAGE(DenseLat3.fuzzy_equals(DenseLat3_mixed,test_precision<data_type>()),std::string("Full Dimension Solve Test 2 for ")+typeid(data_type).name());



    DenseLat3=DenseLat1_lsqr.solve(DenseLat2);
    DenseLat3_mixed=DenseLat1_lsqr.solve(SparseLat2);

    BOOST_CHECK_MESSAGE(DenseLat3.fuzzy_equals(DenseLat3_mixed,test_precision<data_type>()),std::string("LSQR Solve Test 1 for ")+typeid(data_type).name());



}

BOOST_AUTO_TEST_CASE( MixedLatticeSolveTests )
{


    //multwork<double>(3,3,3,3);
    //solvework<double>(5,5,5,5);
    solvework<double>(20,20,20,20);
    solvework<float>(20,20,20,20);




}
