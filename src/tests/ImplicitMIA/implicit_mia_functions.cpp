#include <iostream>
#include <typeinfo>
#define BOOST_TEST_MODULE ImplicitMIAFunctionTests



#include "MIAConfig.h"

#ifdef MIA_USE_HEADER_ONLY_TESTS
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif


#include "ImplicitMIA.h"
#include "DenseMIA.h"
#include "LibMIAUtil.h"
#include "Index.h"

template<class _data_type>
void functions_work(){




    size_t dim1=2;
    size_t dim2=3;
    size_t dim3=4;

    typedef LibMIA::ImplicitMIA<_data_type,3> MIAType;
    typedef LibMIA::ImplicitMIA<_data_type,3,true> MIAType_ref;
    typedef LibMIA::DenseMIA<_data_type,3> DenseMIAType;
    typedef typename LibMIA::internal::index_type<MIAType>::type index_type;
    typedef typename LibMIA::internal::function_type<MIAType>::type function_type;
    typedef typename LibMIA::internal::function_type<MIAType_ref>::type function_type_ref;
    function_type _func=[](index_type idx){
        if(idx%2==0)
            return 1;
        else
            return 2;
    };
    MIAType a(_func,dim1,dim2,dim3);
    bool _correct=true;
    for(auto it=a.data_begin();it<a.data_end();++it){
        if((it-a.data_begin())%2==0){
            if(*it!=1)
                _correct=false;
        }
        else{
            if(*it!=2)
                _correct=false;
        }

    }
    BOOST_CHECK_MESSAGE(_correct,std::string("Basic data iterator test for ImplicitMIA for ")+typeid(_data_type).name());

    DenseMIAType dense_a(dim1,dim2,dim3);
    for(size_t idx=0;idx<dense_a.dimensionality();++idx){
        if(idx%2==0)
            dense_a.atIdx(idx)=1;
        else
            dense_a.atIdx(idx)=2;

    }
    BOOST_CHECK_MESSAGE(a==dense_a,std::string("Basic ImplicitMIA Function Test for ")+typeid(_data_type).name());

    //Now test ImplicitMIAs that refer to another piece of data
    function_type_ref _func_ref=[&dense_a](index_type idx)->_data_type&{

        return dense_a.atIdx(idx);;
    };
    MIAType_ref a_ref(_func_ref,dim1,dim2,dim3);
    BOOST_CHECK_MESSAGE(a_ref==dense_a,std::string("Value test for ImplicitMIA with reference for ")+typeid(_data_type).name());
    _correct=true;
    for(size_t idx=0;idx<dense_a.dimensionality();++idx){
        if(&(dense_a.atIdx(idx))!=&(a_ref.atIdx(idx))){
            _correct=false;
            break;
        }

    }
    BOOST_CHECK_MESSAGE(_correct,std::string("Reference test for ImplicitMIA with reference for ")+typeid(_data_type).name());

    //now test the assign to data_type operation - which implicity tests the iterators of ImplicitMIA
    a_ref.fill(2);
    _correct=true;
    for(size_t idx=0;idx<a_ref.dimensionality();++idx){
        if(a_ref.atIdx(idx)!=2){
            _correct=false;
            break;
        }

    }
    BOOST_CHECK_MESSAGE(_correct,std::string("Basic assignment functionality test to referred data for ImplicitMIA with reference for ")+typeid(_data_type).name());
    _correct=true;
    for(size_t idx=0;idx<dense_a.dimensionality();++idx){
        if(dense_a.atIdx(idx)!=2){
            _correct=false;
            break;
        }

    }

}

BOOST_AUTO_TEST_CASE( ImplicitMIAFunctionTests )
{



    functions_work<double>();

    functions_work<float>();
    functions_work<int>();
    functions_work<long>();






}
