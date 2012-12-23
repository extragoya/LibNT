#ifndef UTIL_H
#define UTIL_H

#include<array>

#include <boost/type_traits.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector.hpp>
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
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/begin_end.hpp>

#include "Index.h"

namespace LibMIA
{
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

template <class T,size_t _order>
class DenseMIA;

template<class _MIA,class m_Seq>
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
        boost::is_same<
        typename boost::numeric::conversion_traits<index_type,arg >::supertype,
        index_type
        >,
        check_mia_dim_args<index_type,Args...>
        >
    {};

//checks to make sure size of variadic arguments is equal to order of MIA
template<class _MIA,typename...Args>
struct check_dims_count : boost::mpl::bool_<sizeof...( Args ) ==internal::order<_MIA>::value>
{ };

template<class _MIA, typename...Args>
struct check_mia_constructor
{

    typedef typename
    boost::mpl::and_<
    check_dims_count<_MIA,Args...>,
    check_mia_dim_args<typename internal::index_type<_MIA>::type,Args...>
    >::type type;

};

template<class T> struct incomplete;



//base case - set to true. This only occurs when LSeq is an empty mpl sequence
template<typename LSeq,typename RSeq,bool empty_LSeq,int expression_type,int recursive_depth>
struct cartesian_product_indices
{
    typedef boost::true_type allowed_recursive;
    typedef boost::mpl::vector_c<int> match_order_sequence;
    typedef boost::mpl::vector_c<int> inter_match_order_sequence;
    typedef boost::mpl::vector_c<int> no_match_order_sequence;


};

//this structure checks how many times the first element of LSeq occurs in RSeq
//Once this is done, it then checks to make sure the frequency follows the expression_type rules
//finally it recurses to check the next element of LSeq
template<typename LSeq,typename RSeq,int expression_type,int recursive_depth>
struct cartesian_product_indices<LSeq,RSeq,false,expression_type,recursive_depth>
{

    //Pull first index type off of left sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<LSeq>::type>::type first_LSeq;
    //count its occurences in the right sequence
    typedef typename boost::mpl::count<RSeq,first_LSeq>::type n;
    //pull the allowed matches for LSeq's first index type
    typedef typename internal::match_rule<
    first_LSeq,expression_type
    >::allowed_matches allowed_matches;
    //see if occurence count is found in the allowed matches
    typedef typename boost::mpl::find_if<
    allowed_matches,
    boost::mpl::equal_to<_,n>
    >::type found_location;
    //store whether the number of occurences is allowed or not
    typedef typename boost::mpl::not_<
    boost::is_same<
    found_location,
    typename boost::mpl::end<allowed_matches>::type
    >
    > allowed;
    //pop the first type off of LSeq
    typedef typename boost::mpl::pop_front<LSeq>::type poppedLSeq;
    //obtain next cartesion_product in recursion depth
    typedef  cartesian_product_indices<poppedLSeq,RSeq,boost::mpl::empty<poppedLSeq>::value,expression_type,recursive_depth+1> next_cartesian_product;


    //now check the remaining LSeq types using recursive calls
    typedef typename boost::mpl::and_<
    allowed,
    typename next_cartesian_product::allowed_recursive
    > allowed_recursive;
    constexpr static bool value= allowed_recursive::value;

    typedef typename boost::mpl::if_<
    boost::mpl::and_<
    boost::mpl::equal_to<
    n,boost::mpl::int_<1>
    >,
    boost::mpl::not_<
    boost::mpl::bool_<first_LSeq::elemval>
    >
    >,
    typename boost::mpl::push_front<
    typename next_cartesian_product::match_order_sequence,
    boost::mpl::int_<recursive_depth>
    >::type,
    typename next_cartesian_product::match_order_sequence
    >::type match_order_sequence;


    typedef typename boost::mpl::if_<
    boost::mpl::and_<
    boost::mpl::equal_to<
    n,boost::mpl::int_<0>
    >,
    boost::mpl::not_<
    boost::mpl::bool_<first_LSeq::elemval>
    >
    >,
    typename boost::mpl::push_front<
    typename next_cartesian_product::no_match_order_sequence,
    boost::mpl::int_<recursive_depth>
    >::type,
    typename next_cartesian_product::no_match_order_sequence
    >::type no_match_order_sequence;

    typedef typename boost::mpl::if_<
    boost::mpl::and_<
    boost::mpl::equal_to<
    n,boost::mpl::int_<1>
    >,
    boost::mpl::bool_<first_LSeq::elemval>
    >,
    typename boost::mpl::push_front<
    typename next_cartesian_product::inter_match_order_sequence,
    boost::mpl::int_<recursive_depth>
    >::type,
    typename next_cartesian_product::inter_match_order_sequence
    >::type inter_match_order_sequence;

};


//base case - set to true. This only occurs when LSeq is an empty mpl sequence
template<typename Seq,bool empty_LSeq,int expression_type>
struct auto_cartesian_product_indices
{
    typedef boost::true_type allowed_recursive;
};

//this structure checks how many times the first element of LSeq occurs in RSeq
//Once this is done, it then checks to make sure the frequency follows the expression_type rules
//finally it recurses to check the next element of LSeq
template<typename Seq,int expression_type>
struct auto_cartesian_product_indices<Seq,false,expression_type>
{

    //Pull first index type off of the sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<Seq>::type>::type first_Seq;
    //pop the first type off of Seq
    typedef typename boost::mpl::pop_front<Seq>::type poppedSeq;
    //count its occurences in the sequence
    typedef typename boost::mpl::count<poppedSeq,first_Seq>::type n;
    //pull the allowed matches for Seq's first index type
    typedef typename internal::auto_match_rule<
    first_Seq,expression_type
    >::allowed_matches allowed_matches;
    //see if occurence count is found in the allowed matches
    typedef typename boost::mpl::find_if<
    allowed_matches,
    boost::mpl::equal_to<_,n>
    >::type found_location;
    //store whether the number of occurences is allowed or not
    typedef typename boost::mpl::not_<
    boost::is_same<
    found_location,
    typename boost::mpl::end<allowed_matches>::type
    >
    > allowed;

    //now check the remaining LSeq types using recursive calls
    typedef typename boost::mpl::and_<
    allowed,
    typename auto_cartesian_product_indices<poppedSeq,boost::mpl::empty<poppedSeq>::value,expression_type>::allowed_recursive
    > allowed_recursive;
    constexpr static bool value= allowed_recursive::value;

};



template<typename LSeq, typename RSeq,bool empty_LSeq>
struct pull_product_indices
{

    typedef boost::mpl::vector<> outer_product_indices;
    typedef boost::mpl::vector<> inter_product_indices;


};

template<typename LSeq, typename RSeq>
struct pull_product_indices<LSeq,RSeq,false>
{

    //Pull first index type off of left sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<LSeq>::type>::type first_LSeq;
    //count its occurences in the right sequence
    typedef typename boost::mpl::count<RSeq,first_LSeq>::type n;


    typedef typename boost::mpl::equal_to<boost::mpl::int_<1>,n>::type matched;

    //pop the first type off of LSeq
    typedef typename boost::mpl::pop_front<LSeq>::type poppedLSeq;
    //obtain next cartesion_product in recursion depth
    typedef  pull_product_indices<poppedLSeq,RSeq,boost::mpl::empty<poppedLSeq>::value> next_pull_product_indices;


    typedef typename boost::mpl::if_<
    boost::mpl::and_<
    matched,
    boost::mpl::bool_<first_LSeq::elemval>
    >,
    typename boost::mpl::push_front<
    typename next_pull_product_indices::inter_product_indices,
    first_LSeq
    >::type,
    typename next_pull_product_indices::inter_product_indices
    >::type inter_product_indices;

    typedef typename boost::mpl::if_<
    boost::mpl::and_<
    boost::mpl::not_<
    matched
    >,
    boost::mpl::not_<
    boost::mpl::bool_<first_LSeq::elemval>
    >
    >,
    typename boost::mpl::push_front<
    typename next_pull_product_indices::outer_product_indices,
    first_LSeq
    >::type,
    typename next_pull_product_indices::outer_product_indices
    >::type outer_product_indices;

};




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


template<typename _MIA, typename...Args>
struct check_mia_constructor;


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

template<class L_MIA, class R_MIA, unsigned int order>
struct MIAProductReturnType
{
};

template<class Derived,class otherDerived, unsigned int order>
struct MIAProductReturnType<DenseMIABase<Derived>,DenseMIABase<otherDerived>,order>
{
    typedef DenseMIA<typename ScalarPromoteType<Derived,otherDerived>::type,order> MIA_return_type;

};


template<class L_MIA, class R_MIA,class l_Seq,class r_Seq>
struct MIAProductUtil
{



    typedef typename internal::pull_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value>::outer_product_indices l_indices;
    typedef typename internal::pull_product_indices<r_Seq,l_Seq,boost::mpl::empty<r_Seq>::value>::outer_product_indices r_indices;
    typedef typename internal::pull_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value>::inter_product_indices inter_product_indices;
    typedef typename boost::mpl::insert_range<l_indices,typename boost::mpl::end<l_indices>::type,r_indices>::type concat;

    typedef typename boost::mpl::insert_range<concat,typename boost::mpl::end<concat>::type,inter_product_indices>::type final_sequence;
    static constexpr int MIA_return_order= boost::mpl::size<final_sequence>::value;
    typedef typename MIAProductReturnType<L_MIA,R_MIA,MIA_return_order>::MIA_return_type MIA_return_type;

};

//check, probably std::array second parameter templated by size_t
template<typename indexType,size_t T>
std::array<indexType,T> ind2sub(indexType i, std::array<indexType,T> dims)
{

    std::array<indexType,T> indices;
    indexType divisor=1;
    for(size_t i=0; i<T; ++i)
    {
        indices[i]=i/divisor;
        indices[i]%=dims[i];
        divisor*=dims[i];
    }

}

template<typename indexType,size_t T>
indexType sub2ind(std::array<indexType,T> indices, std::array<indexType,T> dims)
{


    indexType idx=0;
    indexType multiplier=1;
    for(size_t i=0; i<T; ++i)
    {
        idx+=indices[i]*multiplier;
        multiplier*=dims[i];
    }

}

template<typename accessType, typename accessType2,typename dimType>
typename dimType::value_type sub2ind(accessType indices, accessType2 order, dimType dims)
{



    typename dimType::value_type idx=0;
    typename dimType::value_type multiplier;

    for (size_t i=0; i< indices.size(); ++i)
    {
        multiplier=1;
        for(size_t j=0; j<order[i]; ++j)
        {

            multiplier*=dims[j];
        }
        idx+=indices[i]*multiplier;

    }


}

/*! @} */
/*! @} */

}


#endif // SPARSEUTIL_H
