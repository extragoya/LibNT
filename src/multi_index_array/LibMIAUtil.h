// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.

//This file should be dedicated primarily to forward declare and datatype resolutions

#ifndef LIBMIAUTIL_H
#define LIBMIAUTIL_H


#include <time.h>


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>


#define PARALLEL_TOL 8192
//#define LIBMIA_CHECK_DIMS 1



namespace LibMIA
{

inline boost::mt19937& LibMIA_gen() {
   static boost::mt19937 gen(time(0));
   return gen;
}

/** \addtogroup util Utilities
 *  @{
*/

/** \addtogroup lattice_util Lattice Utilities
 *  @{
 */



enum SolveInfo{
    FullyRanked,
    LeastSquares,
    NoInfo,
    RankDeficient
};

const bool ColumnMajor=true;
const bool RowMajor=false;

//workaround for MSVC compiler which doesn't yet support constexpr
#ifdef NO_CONST_EXPR_LIBMIA

//comparison tolerances for non zeros after operations like solving
template <class T> struct Tolerance
{
	const static int tolerance;
};
template <class T>
const int __declspec(selectany)  Tolerance<T>::tolerance = 0;
template <> struct Tolerance<float>
{
	const static float tolerance; 
	//static const float tolerance=5.96e-08;

};
const float __declspec(selectany)  Tolerance<float>::tolerance = 1e-3f;
template <> struct Tolerance<double>
{
	const static double tolerance; 
	//static const double tolerance=1.11e-16;

};
const double __declspec(selectany)  Tolerance<double>::tolerance = 1e-6;
//tolerances for how close to zero a nonzero can get to be included in a sparse data structure
template <class T> struct SparseTolerance
{
	const static int tolerance;
};
template <class T>
const int __declspec(selectany)  SparseTolerance<T>::tolerance = 0;
template <> struct SparseTolerance<float>
{
	const static float tolerance;
	//static const float tolerance=5.96e-08;

};
const float __declspec(selectany)  SparseTolerance<float>::tolerance = 5.96e-08;
template <> struct SparseTolerance<double>
{
	const static double tolerance; 
	//static const double tolerance=1.11e-16;

};
const double __declspec(selectany)  SparseTolerance<double>::tolerance = 1.11e-16;
#else


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

#endif


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

template <class T,size_t _order,bool isRef=false>
class ImplicitMIA;

template<class _MIA,class m_Seq,bool ownership=true,size_t inter_product_size=0>
class MIA_Atom;

template<size_t ID,int ElemWise>
struct ProdInd;

template<typename index_type=long long>
struct Range;

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

template<class T, size_t _order,bool isRef>
struct is_MIA<ImplicitMIA<T,_order,isRef > >: public boost::true_type {};

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

template<class T, size_t _order,bool isRef>
struct is_DenseMIA<ImplicitMIA<T,_order,isRef > >: public boost::true_type {};

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

template<class T>
struct is_ProdInd: public boost::false_type {};

template<size_t i,int elemval>
struct is_ProdInd<ProdInd<i,elemval>>: public boost::true_type {};

template<class T>
struct is_Range: public boost::false_type {};

template<class index_type>
struct is_Range<Range<index_type>>: public boost::true_type {};

//see below
template<class T,class Enable = void>
struct const_argument_qualifier
{
    typedef const T type;
};

//SparseMIAs and Lattices could be modified (by a re-sort) even when the operation doesn't alter the underlying data. So in C++, even if the underlying data in 'const', the datatype will not be
template<class T>
struct const_argument_qualifier<T,typename boost::enable_if<is_SparseMIA<T> >::type>
{
    typedef T type;
};

template<class T>
struct const_argument_qualifier<T,typename boost::enable_if<is_SparseLattice<T> >::type>
{
    typedef T type;
};


template<class T> struct incomplete;









//should be undefined
template<typename...Args>
struct check_mia_index_args;

//base case, inherit from true_type
template<>
struct check_mia_index_args<>:public boost::true_type {};

//checks to ensure dimension arguments are all convertable to index_type. Uses recursion to allow
//check to happen regardless of number of arguments. Checking is controlled through inheritance
template<typename arg,typename...Args>
struct check_mia_index_args<arg,Args...> :
        boost::mpl::and_<
            internal::is_ProdInd<arg>,
            check_mia_index_args<Args...>
        >
    {};


//should be undefined
template<class index_type,typename...Args>
struct check_mia_dim_args;


//base case, inherit from true_type
template<class index_type>
struct check_mia_dim_args<index_type>:public boost::true_type {};

//checks to ensure dimension arguments are all convertable to index_type. Uses recursion to allow
//check to happen regardless of number of arguments. Checking is controlled through inheritance
template<class index_type,typename arg,typename...Args>
struct check_mia_dim_args<index_type,arg,Args...> :
        boost::mpl::and_<
        std::is_integral<arg>,
        check_mia_dim_args<index_type,Args...>
        >
    {};

//checks to make sure size of variadic arguments is equal to order of MIA
template<class _MIA,typename...Args>
struct check_dims_count : boost::mpl::bool_<sizeof...( Args ) ==internal::order<_MIA>::value>
{ };

//checks to make sure size of variadic arguments is equal to order of MIA
template<class _MIA>
struct check_order : boost::mpl::bool_<(internal::order<_MIA>::value >0)>
{ };

template<typename _MIA, typename...Args>
struct check_mia_constructor;

template<class _MIA, typename...Args>
struct check_mia_constructor
{

    typedef typename
    boost::mpl::and_<
        boost::mpl::and_<
            check_dims_count<_MIA,Args...>,
            check_mia_dim_args<typename internal::index_type<_MIA>::type,Args...>
        >,
        check_order<_MIA>
    >::type type;

};

template<typename _MIA, typename...Args>
struct check_mia_indexing;

template<class _MIA, typename...Args>
struct check_mia_indexing
{

    typedef typename
    boost::mpl::and_<
        check_dims_count<_MIA,Args...>,
        check_mia_index_args<Args...>
    >::type type;

};

//should be undefined
template<typename...Ranges>
struct check_mia_range_args;

//base case, inherit from true_type
template<>
struct check_mia_range_args<>:public boost::true_type {};

//checks to ensure that ranges are either Range type or integral type. Uses recursion to allow
//check to happen regardless of number of arguments. Checking is controlled through inheritance
template<typename range,typename...Ranges>
struct check_mia_range_args<range,Ranges...> :
        boost::mpl::and_<
            boost::mpl::or_<
                internal::is_Range<range>,
                boost::is_integral<range>
            >,
            check_mia_range_args<Ranges...>
        >
    {};




//checks that the number of Range arguments equals the MIA order and that the arguments are the appropriate data type
template<class _MIA, typename...Ranges>
struct check_ranges
{

    typedef typename
    boost::mpl::and_<
        check_dims_count<_MIA,Ranges...>,
        check_mia_range_args<Ranges...>
    >::type type;

};

template<typename...Ranges>
struct get_range_count{
    typedef typename Indicial_Sequence<Ranges...>::sequence mpl_ranges;
    typedef typename
        boost::mpl::count_if<
            mpl_ranges,
            internal::is_Range<_>
        >::type range_count;

};




template<class index_type1,class index_type2>
struct check_index_compatibility:
    std::is_integral<index_type2>
{};









}

//determines the super datatype in Lhs and Rhs
template<typename Lhs, typename Rhs>
struct ScalarPromoteType
{
    typedef typename boost::numeric::conversion_traits<typename internal::data_type<typename std::remove_const<Lhs>::type>::type,typename internal::data_type<typename std::remove_const<Rhs>::type >::type>::supertype type;

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
    boost::is_same<
    typename internal::data_type<Lhs>::type,
    typename internal::data_type<Rhs>::type
    >,
    SparseLattice<typename internal::data_type<Lhs>::type>
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
            internal::is_SparseMIA<typename std::remove_const<L_MIA>::type>,
            internal::is_SparseMIA<typename std::remove_const<R_MIA>::type>
        >
    >::type
>
{
    typedef SparseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

template<class L_MIA,class R_MIA, size_t order
>
struct MIAProductReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::and_<
            internal::is_DenseMIA<typename std::remove_const<L_MIA>::type>,
            internal::is_DenseMIA<typename std::remove_const<R_MIA>::type>
        >
    >::type
>
{
    typedef DenseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

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
            internal::is_DenseMIA<typename std::remove_const<L_MIA>::type>,
            internal::is_DenseMIA<typename std::remove_const<R_MIA>::type >
        >
    >::type
>
{
    typedef DenseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

template<class L_MIA,class R_MIA, size_t order
>
struct MIANoLatticeProductReturnType<L_MIA,R_MIA,order,
    typename boost::enable_if<
        boost::mpl::or_< //when no lattice mapping is required, dense * sparse is always a sparse
            internal::is_SparseMIA<typename std::remove_const<L_MIA>::type>,
            internal::is_SparseMIA<typename std::remove_const<R_MIA>::type>
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
            internal::is_MIA<typename std::remove_const<L_MIA>::type>,
            internal::is_MIA<typename std::remove_const<R_MIA>::type>
        >
    >::type
>
{
    typedef DenseMIA<typename ScalarPromoteType<L_MIA,R_MIA>::type,order> type;

};

//the datatype after application of a nonlinear function (copied).
template<typename MIA,class Enable=void>
struct MIANonlinearFuncType{};

template<typename MIA>
struct MIANonlinearFuncType<MIA,typename boost::enable_if<internal::is_DenseMIA<typename std::remove_const<MIA>::type>>::type>{
    typedef typename std::remove_const<MIA>::type MIA_type;
    typedef DenseMIA<typename internal::data_type<MIA_type>::type,internal::order<MIA_type>::value> type;
};

template<typename MIA>
struct MIANonlinearFuncType<MIA,typename boost::enable_if<internal::is_SparseMIA<typename std::remove_const<MIA>::type>>::type>{
    typedef typename std::remove_const<MIA>::type MIA_type;
    typedef SparseMIA<typename internal::data_type<MIA_type>::type,internal::order<MIA_type>::value> type;

};


template<typename MIA,size_t remaining_indices_size,class Enable = void>
struct MIAUnaryType
{

};


template<typename MIA,size_t remaining_indices_size>
struct MIAUnaryType<MIA,remaining_indices_size,typename boost::enable_if<internal::is_DenseMIA<typename std::remove_const<MIA>::type>>::type>
{
    typedef DenseMIA<typename internal::data_type<typename std::remove_const<MIA>::type>::type,remaining_indices_size> type;
};

template<typename MIA,size_t remaining_indices_size>
struct MIAUnaryType<MIA,remaining_indices_size,typename boost::enable_if<internal::is_SparseMIA<typename std::remove_const<MIA>::type>>::type>
{
    typedef SparseMIA<typename internal::data_type<typename std::remove_const<MIA>::type>::type,remaining_indices_size> type;
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
                internal::is_DenseMIA<typename std::remove_const<Lhs>::type>,
                internal::is_DenseMIA<typename std::remove_const<Rhs>::type>
            >
        >::type
    >
{

    typedef DenseMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<typename std::remove_const<Lhs>::type>::value> type;
};


template<typename Lhs, typename Rhs>
struct MIAMergeReturnType<Lhs,Rhs,
    typename
        boost::enable_if<
            boost::mpl::or_<
                boost::mpl::and_<
                    internal::is_DenseMIA<typename std::remove_const<Lhs>::type>,
                    internal::is_SparseMIA<typename std::remove_const<Rhs>::type>
                >,
                boost::mpl::and_<
                    internal::is_SparseMIA<typename std::remove_const<Lhs>::type>,
                    internal::is_DenseMIA<typename std::remove_const<Rhs>::type>
                >
            >
        >::type
    >
{

    typedef DenseMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<typename std::remove_const<Lhs>::type>::value> type;
};

template<typename Lhs, typename Rhs>
struct MIAMergeReturnType<Lhs,Rhs,
    typename
        boost::enable_if<
            boost::mpl::and_<
                internal::is_SparseMIA<typename std::remove_const<Lhs>::type>,
                internal::is_SparseMIA<typename std::remove_const<Rhs>::type>
            >
        >::type
    >
{

    typedef SparseMIA<typename ScalarPromoteType<Lhs,Rhs>::type,internal::order<typename std::remove_const<Lhs>::type>::value> type;
};






/*! @} */
/*! @} */

}


#endif // LIBMIAUTIL_H

