// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef UTIL_H
#define UTIL_H

#include<array>
#include <time.h>

#include <boost/type_traits.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/print.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/count.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/print.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/tuple/tuple.hpp>

#include "PermuteIterator.h"
#include "Index.h"

namespace LibMIA
{
boost::random::mt19937 gen(time(0));
/** \addtogroup util Utilities
 *  @{
*/

/** \addtogroup lattice_util Lattice Utilities
 *  @{
 */

template<class T>
struct MIAprint
{
    void operator() (T i)
    {
        std::cout << " " << i;
    }
} ;

template<class T>
struct select_first
{
    T& operator()(T&left, T& right){
        return left;
    }
};



enum SolveInfo{
    FullyRanked,
    LeastSquares,
    NoInfo,
    RankDeficient
};

const bool ColumnMajor=true;
const bool RowMajor=false;


//comparison tolerances for non zeros after operations like solving
template <class T> struct Tolerance
{
    constexpr static int tolerance=0;
};
template <> struct Tolerance<float>
{
    constexpr static float tolerance=1e-3f;
    //static const float tolerance=5.96e-08;

};
template <> struct Tolerance<double>
{
    constexpr static double tolerance=1e-6;
    //static const double tolerance=1.11e-16;

};

//tolerances for how close to zero a nonzero can get to be included in a sparse data structure
template <class T> struct SparseTolerance
{
    constexpr static int tolerance=0;
};
template <> struct SparseTolerance<float>
{
    constexpr static float tolerance=5.96e-08;
    //static const float tolerance=5.96e-08;

};
template <> struct SparseTolerance<double>
{
    constexpr static double tolerance=1.11e-16;
    //static const double tolerance=1.11e-16;

};



inline long double log2(const long double x){
    return  std::log(x) * M_LOG2E;
}


template <class Derived>
class Lattice;

template <class Derived>
class SparseLatticeBase;

template <class T>
class SparseLattice;

template <class T>
class MappedSparseLattice;


template <class Derived>
class DenseLatticeBase;

template <class T>
class DenseLattice;

template <class T>
class MappedDenseLattice;

template <class Derived>
class MIA;

template <class Derived>
class DenseMIABase;

template <class Derived>
class SparseMIABase;

template <class T,size_t _order>
class DenseMIA;

template <class T,size_t _order>
class SparseMIA;

template <class T,size_t _order>
class ImplicitMIA;

template<class _MIA,class m_Seq,bool ownership=true,size_t inter_product_size=0>
class MIA_Atom;



namespace internal
{

using boost::mpl::_;

template<class Derived>
struct data_type;

template<class Derived>
struct data_type_ref;

template<class Derived>
struct const_data_type_ref;

template<class Derived>
struct order;

template<class Derived>
struct index_type;

template<class Derived>
struct Data;

template<typename Derived>
struct Indices;

template<typename Derived>
struct full_iterator_tuple;

template<typename Derived>
struct const_full_iterator_tuple;

template<typename Derived>
struct full_tuple;

template<typename Derived>
struct const_full_tuple;

template<typename Derived>
struct storage_iterator;

template<typename Derived>
struct const_storage_iterator;

template<typename Derived>
struct index_iterator;

template<typename Derived>
struct const_index_iterator;

template<typename Derived>
struct data_iterator;

template<typename Derived>
struct const_data_iterator;

//only use for implicit MIAs
template<typename Derived>
struct function_type;

template<typename Derived>
struct FinalDerived;

template<typename... Ts>
struct Indicial_Sequence;

template<class T>
struct is_MIA: public boost::false_type {};

template<class Derived>
struct is_MIA<MIA<Derived > >: public boost::true_type {};

template<class Derived>
struct is_MIA<DenseMIABase<Derived > >: public boost::true_type {};

template<class T, size_t _order>
struct is_MIA<DenseMIA<T,_order > >: public boost::true_type {};

template<class T, size_t _order>
struct is_MIA<ImplicitMIA<T,_order > >: public boost::true_type {};

template<class T, size_t _order>
struct is_MIA<SparseMIA<T,_order > >: public boost::true_type {};

template<class Derived>
struct is_MIA<SparseMIABase<Derived> >: public boost::true_type {};

template<class T>
struct is_DenseMIA: public boost::false_type {};

template<class T, size_t _order>
struct is_DenseMIA<DenseMIA<T,_order > >: public boost::true_type {};

template<class Derived>
struct is_DenseMIA<DenseMIABase<Derived> >:public boost::true_type {};

template<class T, size_t _order>
struct is_DenseMIA<ImplicitMIA<T,_order > >: public boost::true_type {};

template<class Derived>
struct is_DenseMIA<MIA<Derived> >:public is_DenseMIA<Derived> {};

template<class T>
struct is_SparseMIA: public boost::false_type {};

template<class T, size_t _order>
struct is_SparseMIA<SparseMIA<T,_order > >: public boost::true_type {};

template<class Derived>
struct is_SparseMIA<SparseMIABase<Derived> >: public boost::true_type {};


template<class T>
struct is_Lattice: public boost::false_type {};

template<class T>
struct is_SparseLattice: public boost::false_type {};

template<class T>
struct is_DenseLattice: public boost::false_type {};

template<class Derived>
struct is_Lattice<Lattice<Derived > >: public boost::true_type {};

template<class Derived>
struct is_SparseLattice<SparseLatticeBase<Derived > >: public boost::true_type {};

template<class Derived>
struct is_DenseLattice<DenseLatticeBase<Derived > >: public boost::true_type {};



template<class T> struct incomplete;

//! Converts a scalar value to data_type
/*!
    \tparam from_data_type the data_type you are converting from
*/
template<class data_type,class from_data_type,typename boost::enable_if< boost::is_pod< from_data_type >, int >::type = 0>
data_type convert(const from_data_type from){
    using namespace boost::numeric;
    typedef boost::numeric::converter<data_type,boost::uniform_real<>::result_type> to_mdata_type;
    return to_mdata_type::convert(from);
}




//must be boost::tuples of iterators. Assumes a's container is sized to be a.size+b.size
template<class ADataIt, class AIndexIt, class BDataIt, class BIndexIt,class Op>
ADataIt merge_sparse_storage_containers(ADataIt  a_data_begin,ADataIt  a_data_end,AIndexIt  a_index_begin,AIndexIt  a_index_end,BDataIt  b_data_begin,BDataIt  b_data_end,BIndexIt  b_index_begin,BIndexIt  b_index_end,Op op)
{
    using namespace boost::numeric;
    typedef typename ADataIt::value_type a_data_type;
    typedef typename BDataIt::value_type b_data_type;


    typedef converter<a_data_type,b_data_type,conversion_traits<a_data_type,b_data_type>,def_overflow_handler,RoundEven<b_data_type>> to_mdata_type;
    ADataIt a_actual_data_end=a_data_end;
    AIndexIt a_actual_index_end=a_index_end;
    ADataIt a_cur_data_it=a_data_begin;
    AIndexIt a_cur_index_it=a_index_begin;
    while(a_cur_data_it<a_data_end && b_data_begin<b_data_end){
        if (*a_cur_index_it<*b_index_begin){
            a_cur_index_it++;
            a_cur_data_it++;
        }
        else if  (*b_index_begin<*a_cur_index_it){
            *a_actual_data_end++=*b_data_begin++;
            *a_actual_index_end++=*b_index_begin++;


        }
        else{
            *a_cur_data_it=op(*a_cur_data_it,to_mdata_type::convert(*b_data_begin++));
            a_cur_data_it++;
            a_cur_index_it++;
            b_index_begin++;
        }

    }
    if (a_cur_data_it==a_data_end){
        while (b_data_begin<b_data_end){
            *a_actual_data_end++=*b_data_begin++;
            *a_actual_index_end++=*b_index_begin++;

        }
    }
    std::cout << "Index\t Data in scan merge" << std::endl;
    for(auto i=a_data_begin,j=a_index_begin;i<a_actual_data_end;++i,++j)
        std::cout << *j << "\t " << *i << std::endl;

    std::cout << std::endl;

    std::cout << " diff " << a_index_end-a_index_begin << std::endl;
    std::inplace_merge(make_sort_permute_iter(a_index_begin,a_data_begin),
                    make_sort_permute_iter(a_index_end,a_data_end),
              make_sort_permute_iter(a_actual_index_end,a_actual_data_end),
                    sort_permute_iter_compare<AIndexIt,ADataIt>());

//    std::inplace_merge(a_data_begin,a_data_end,a_actual_data_end,[&](const typename ADataIt::value_type& lhs, const typename ADataIt::value_type& rhs)
//    {
//        return *(a_index_begin+(&lhs-&(*a_data_begin))) <*(a_index_begin+(&rhs-&(*a_data_begin)));
//    });
//    std::inplace_merge(a_index_begin,a_index_end,a_actual_index_end);
    std::cout << " diff " << a_index_end-a_index_begin << std::endl;
    std::cout << "Index\t Data in AFTER scan merge" << std::endl;
    for(auto i=a_data_begin,j=a_index_begin;i<a_actual_data_end;++i,++j)
        std::cout << *j << "\t " << *i << std::endl;

    std::cout << std::endl;
    return a_actual_data_end;

}


//must be boost::tuples of iterators. Assumes a's container is sized to be a.size+b.size
template<class AStorageItType, class BStorageItType, class Op>
AStorageItType merge_sparse_storage_containers(AStorageItType  a_begin,AStorageItType  a_end,BStorageItType  b_begin,BStorageItType  b_end,Op op)
{
    using namespace boost::numeric;
    typedef typename boost::remove_reference<typename BStorageItType::value_type::first_type>::type b_data_type;
    typedef typename boost::remove_reference<typename AStorageItType::value_type::first_type>::type a_data_type;

    typedef converter<a_data_type,b_data_type,conversion_traits<a_data_type,b_data_type>,def_overflow_handler,RoundEven<b_data_type>> to_mdata_type;
    AStorageItType a_actual_end=a_end;
    AStorageItType a_actual_begin=a_begin;
    while(a_begin<a_end && b_begin<b_end){
        if (std::get<1>(*a_begin)<std::get<1>(*b_begin)){
            a_begin++;
        }
        else if  (std::get<1>(*b_begin)<std::get<1>(*a_begin)){
            std::get<0>(*a_actual_end)=op(a_data_type(0),to_mdata_type::convert(std::get<0>(*b_begin)));
            std::get<1>(*a_actual_end++)=std::get<1>(*b_begin++);


        }
        else{
            std::get<0>(*a_begin)=op(std::get<0>(*a_begin),to_mdata_type::convert(std::get<0>(*b_begin)));
            a_begin++;
            b_begin++;
        }

    }
    if (a_begin==a_end){
        while (b_begin<b_end){
            std::get<0>(*a_actual_end)=op(a_data_type(0),to_mdata_type::convert(std::get<0>(*b_begin)));
            std::get<1>(*a_actual_end++)=std::get<1>(*b_begin++);
        }
    }

    std::inplace_merge(a_actual_begin,a_end,a_actual_end,[](const typename AStorageItType::value_type& lhs, const typename AStorageItType::value_type& rhs)
    {
        return std::get<1>(lhs)<std::get<1>(rhs);
    });


    return a_actual_end;

}

//assumes C, A, and B are of the same dimensions and in the same sort order (and A and B are sorted)
template<class ADerived, class BDerived, class c_data_type,size_t c_order,class Op>
void outside_merge_sparse_storage_containers(SparseMIA<c_data_type,c_order> & C , const SparseMIABase<ADerived> & A,const SparseMIABase<BDerived> & B,Op op)
{
    using namespace boost::numeric;
    typedef typename ADerived::data_type a_data_type;

    C.clear();
    C.reserve(A.size()+B.size());
    auto a_begin=A.index_begin();
    auto b_begin=B.index_begin();
    auto a_end=A.index_end();
    auto b_end=B.index_end();
    while(a_begin<a_end && b_begin<b_end){
        if (*a_begin<*b_begin){
            C.push_back(C.convert(A.data_at(a_begin-A.index_begin())),*a_begin);
            a_begin++;
        }
        else if  (*b_begin<*a_begin){
            C.push_back(C.convert(op(a_data_type(0),B.data_at(b_begin-B.index_begin()))),*b_begin);
            b_begin++;
        }
        else{
            C.push_back(C.convert(op(A.data_at(a_begin-A.index_begin()),B.data_at(b_begin-B.index_begin()))),*a_begin);
            a_begin++;
            b_begin++;
        }

    }
    if (a_begin==a_end){
        while (b_begin<b_end){
            C.push_back(C.convert(op(a_data_type(0),B.data_at(b_begin-B.index_begin()))),*b_begin);
            b_begin++;
        }
    }
    else{
        while (a_begin<a_end){
            C.push_back(C.convert(A.data_at(a_begin-A.index_begin())),*a_begin);
            a_begin++;
        }
    }




}





}

//determines the super datatype in Lhs and Rhs
template<typename Lhs, typename Rhs>
struct ScalarPromoteType
{
    typedef typename boost::numeric::conversion_traits<typename internal::data_type<Lhs>::type,typename internal::data_type<Rhs>::type>::supertype type;

};

//determines the super index type in Lhs and Rhs
template<typename Lhs, typename Rhs>
struct IndexPromoteType
{
    typedef typename boost::numeric::conversion_traits<typename internal::index_type<Lhs>::type,typename internal::index_type<Rhs>::type>::supertype type;

};





//sparse sparse lattice product only enabled if operands share both data and index types
template<typename Lhs, typename Rhs>
struct SparseProductReturnType
{

    typedef typename
    boost::enable_if<
    boost::mpl::and_<
    boost::is_same<
    typename internal::index_type<Lhs>::type,
    typename internal::index_type<Rhs>::type
    >,
    boost::is_same<
    typename internal::data_type<Lhs>::type,
    typename internal::data_type<Rhs>::type
    >
    >
    ,  SparseLattice<typename internal::data_type<Lhs>::type>
    >::type type;


};

//sparse sparse lattice product only enabled if  LHS data type is floating point
template<typename Lhs, typename Rhs>
struct SparseSolveReturnType
{

    typedef
    typename boost::enable_if<
    boost::is_floating_point<
    typename internal::data_type<Lhs>::type
    >, //floating point datatypes
    DenseLattice<typename internal::data_type<Lhs>::type> //return type
    >::type type;
};

template<typename Lhs, typename Rhs>
struct SparseMergeReturnType
{

    typedef SparseLattice<typename ScalarPromoteType<Lhs,Rhs>::type> type;
};

template<typename Lhs, typename Rhs>
struct DenseMergeReturnType
{

    typedef DenseLattice<typename ScalarPromoteType<Lhs,Rhs>::type> type;
};

template<typename Lhs, typename Rhs>
struct DenseProductReturnType
{

    typedef
    typename boost::enable_if<
    boost::is_same<
    typename internal::data_type<Lhs>::type,
    typename internal::data_type<Rhs>::type
    >
    ,  DenseLattice<typename ScalarPromoteType<Lhs,Rhs>::type >
    >::type type;

};

template<typename Lhs, typename Rhs>
struct DenseSolveReturnType
{

    typedef
    typename boost::enable_if<
        boost::is_floating_point<
            typename internal::data_type<Lhs>::type
        > //floating point datatypes
        , typename DenseProductReturnType<Lhs,Rhs>::type //return type
    >::type type;
};

template<class L_MIA, class R_MIA, size_t order,class Enable = void>
struct MIAProductReturnType
{
};



template<class L_MIA,class R_MIA, size_t order
>
struct MIAProductReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::or_< //when a lattice mapping is required, dense * sparse is always a dense
            internal::is_DenseMIA<L_MIA>,
            internal::is_DenseMIA<R_MIA>
        >
    >::type
>
{
    typedef DenseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

template<class L_MIA,class R_MIA, size_t order
>
struct MIAProductReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::and_<
            internal::is_SparseMIA<L_MIA>,
            internal::is_SparseMIA<R_MIA>
        >
    >::type
>
{
    typedef SparseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};


template<class L_MIA, class R_MIA, size_t order,class Enable = void>
struct MIANoLatticeProductReturnType
{
};



template<class L_MIA,class R_MIA, size_t order
>
struct MIANoLatticeProductReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::and_<
            internal::is_DenseMIA<L_MIA>,
            internal::is_DenseMIA<R_MIA>
        >
    >::type
>
{
    typedef ImplicitMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

template<class L_MIA,class R_MIA, size_t order
>
struct MIANoLatticeProductReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::or_< //when no lattice mapping is required, dense * sparse is always a sparse
            internal::is_SparseMIA<L_MIA>,
            internal::is_SparseMIA<R_MIA>
        >
    >::type
>
{
    typedef SparseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

template<class L_MIA, class R_MIA, size_t order,class Enable = void>
struct MIASolveReturnType
{
};
template<class L_MIA,class R_MIA, size_t order>
struct MIASolveReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::and_<
            internal::is_MIA<L_MIA>,
            internal::is_MIA<R_MIA>
        >
    >::type
>
{
    typedef DenseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

template<typename Lhs, typename Rhs,class Enable = void>
struct MIAMergeReturnType
{
};

template<typename Lhs, typename Rhs>
struct MIAMergeReturnType<Lhs,Rhs,
    typename
        boost::enable_if<
            boost::mpl::and_<
                internal::is_DenseMIA<Lhs>,
                internal::is_DenseMIA<Rhs>
            >
        >::type
    >
{

    typedef ImplicitMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<Lhs>::value> type;
};


template<typename Lhs, typename Rhs>
struct MIAMergeReturnType<Lhs,Rhs,
    typename
        boost::enable_if<
            boost::mpl::or_<
                boost::mpl::and_<
                    internal::is_DenseMIA<Lhs>,
                    internal::is_SparseMIA<Rhs>
                >,
                boost::mpl::and_<
                    internal::is_SparseMIA<Lhs>,
                    internal::is_DenseMIA<Rhs>
                >
            >
        >::type
    >
{

    typedef DenseMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<Lhs>::value> type;
};

template<typename Lhs, typename Rhs>
struct MIAMergeReturnType<Lhs,Rhs,
    typename
        boost::enable_if<
            boost::mpl::and_<
                internal::is_SparseMIA<Lhs>,
                internal::is_SparseMIA<Rhs>
            >
        >::type
    >
{

    typedef SparseMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<Lhs>::value> type;
};

template<class array_type>
void print_array(const array_type & _array, const std::string &header){
    std::cout << header;
    for(auto & _i:_array){
        std::cout << " " << _i;
    }
    std::cout << std::endl;

}

template<class T1, class T2,size_t _size>
bool compare_arrays(const std::array<T1,_size> & array1, const std::array<T2,_size> & array2){
    typedef boost::numeric::converter<T1,T2> to_mdata_type;
    for(size_t i=0;i<_size;++i)
        if (array1[i]!=to_mdata_type::convert(array2[i]))
            return false;

    return true;


}

template<class data_type>
struct array_converter
{

    template<class other_data_type,size_t _size>
    static std::array<data_type,_size> convert(const std::array<other_data_type,_size> & _from)
    {
        typedef boost::numeric::converter<data_type,other_data_type> to_mdata_type;
        std::array<data_type,_size> ret;
        for(size_t i=0;i<_size;++i)
            ret[i]=to_mdata_type::convert(_from[i]);

        return ret;
    }

    template<size_t _size>
    static std::array<data_type,_size> convert(std::array<data_type,_size> & _from){
        return _from;
    }
};

//!prec must be positive
template<typename T, typename T2,typename T3>
inline bool isEqualFuzzy(T a, T2 b, T3 prec = Tolerance<T>::tolerance)
{
  if(std::abs(a) < 1 || std::abs(b) < 1)
    return std::abs(a-b)<=prec;
  else{
    return std::abs(a - b) <= std::min(std::abs(a), std::abs(b)) * prec;
  }
}




/*! @} */
/*! @} */

}


#endif // SPARSEUTIL_H

