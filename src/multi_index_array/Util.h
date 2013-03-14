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

const bool ColumnMajor=true;
const bool RowMajor=false;

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

template<class _MIA,class m_Seq,size_t inter_product_size=0>
struct MIA_Atom;


namespace internal
{

using boost::mpl::_;

template<class Derived>
struct data_type;

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
struct is_MIA<SparseMIA<T,_order > >: public boost::true_type {};

template<class Derived>
struct is_MIA<SparseMIABase<Derived> >: public boost::true_type {};

template<class T>
struct is_DenseMIA: public boost::false_type {};

template<class T, size_t _order>
struct is_DenseMIA<DenseMIA<T,_order > >: public boost::true_type {};

template<class Derived>
struct is_DenseMIA<DenseMIABase<Derived> >:public boost::true_type {};

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



//must be boost::tuples of iterators. Assumes a's container is sized to be a.size+b.size
template<class AStorageItType, class BStorageItType, class Op>
AStorageItType merge_sparse_storage_containers(AStorageItType  a_begin,AStorageItType  a_end,BStorageItType  b_begin,BStorageItType  b_end,Op op)
{
    using namespace boost::numeric;
    typedef typename boost::remove_reference<typename boost::tuples::element<0,typename BStorageItType::value_type>::type>::type b_data_type;
    typedef typename boost::remove_reference<typename boost::tuples::element<0,typename AStorageItType::value_type>::type>::type a_data_type;

    typedef converter<a_data_type,b_data_type,conversion_traits<a_data_type,b_data_type>,def_overflow_handler,RoundEven<b_data_type>> to_mdata_type;
    AStorageItType a_actual_end=a_end;
    while(a_begin<a_end && b_begin<b_end){
        if (boost::get<1>(*a_begin)<boost::get<1>(*b_begin)){
            a_begin++;
        }
        else if  (boost::get<1>(*b_begin)<boost::get<1>(*a_begin)){
            *a_actual_end=*b_begin;
            a_actual_end++;
            b_begin++;

        }
        else{
            boost::get<0>(*a_begin)=op(boost::get<0>(*a_begin),to_mdata_type::convert(boost::get<0>(*b_begin)));
            a_begin++;
            b_begin++;
        }

    }
    if (a_begin==a_end){
        while (b_begin<b_end){
            *a_actual_end=*b_begin;
            b_begin++;
        }
    }


    std::inplace_merge(a_begin,a_end,a_actual_end,[](const typename AStorageItType::value_type& lhs, const typename AStorageItType::value_type& rhs)
    //std::inplace_merge(a_begin,a_end,a_actual_end,[](const decltype(*a_begin)& lhs, const decltype(*a_begin)& rhs)
    {
        return boost::get<1>(lhs)<boost::get<1>(rhs);
    });
    return a_actual_end;

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
    ,  Lhs
    >::type type;


};

//sparse sparse lattice product only enabled if operands share index types and LHS data type is floating point
template<typename Lhs, typename Rhs>
struct SparseSolveReturnType
{

    typedef
    typename boost::enable_if<
    boost::mpl::and_<
    boost::is_floating_point<
    typename internal::data_type<Lhs>::type
    >, //floating point datatypes
    boost::is_same<
    typename internal::index_type<Lhs>::type,
    typename internal::index_type<Rhs>::type
    >
    >
    , DenseLattice<typename internal::data_type<Lhs>::type> //return type
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


//only enable when Derived and otherDerived are MIAs - needs to be extended to detect for dense vs sparse mias
template<class L_MIA,class R_MIA, size_t order
>
struct MIAProductReturnType<L_MIA,R_MIA,order,
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
                internal::is_MIA<Lhs>,
                internal::is_MIA<Rhs>
            >
        >::type
    >
{

    typedef DenseMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<Lhs>::value> type;
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


/*! @} */
/*! @} */

}


#endif // SPARSEUTIL_H

