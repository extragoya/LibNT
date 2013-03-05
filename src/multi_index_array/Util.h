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


const bool ColumnMajor=true;
const bool RowMajor=false;

template <class T> struct SparseTolerance
{

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
template <> struct SparseTolerance<int>
{
    constexpr static int tolerance=1;
    //static const int tolerance=1;
};
template <> struct SparseTolerance<long>
{
    constexpr static long tolerance=1;
    //static const long tolerance=1;
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

/*! @} */
/*! @} */

}


#endif // SPARSEUTIL_H

