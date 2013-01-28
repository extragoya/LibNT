// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef INDEXUTIL_H
#define INDEXUTIL_H

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

/** \addtogroup index_util Index Utilities
 *  @{
 */


namespace internal
{

template<class MIA_type,class array_type1, class array_type2>
void collect_dimensions(const MIA_type& _mia, const array_type1& _sequence_order,array_type2& _dims,size_t& curIdx)
{

    for(auto _order: _sequence_order)
    {

        _dims[curIdx++]=_mia.dim((size_t)_order);
    }

}

template<class MIA_type,class array_type1, class array_type2>
void collect_dimensions(const MIA_type& _mia, const array_type1& _sequence_order,array_type2& _dims)
{

    size_t curIdx=0;
    collect_dimensions(_mia,_sequence_order,_dims,curIdx);
}

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

template<class index_type1,class index_type2>
struct check_index_compatibility:
    boost::is_same<
        typename boost::numeric::conversion_traits<index_type1,index_type2 >::supertype,
        index_type1
    >
{};

}

//check, probably std::array second parameter templated by size_t
template<typename indexType,size_t T>
std::array<indexType,T> ind2sub(indexType idx, const std::array<indexType,T> & dims)
{

    std::array<indexType,T> indices;

    for(size_t i=0; i<T; ++i)
    {

        indices[i]=idx%dims[i];
        idx/=dims[i];

    }

    return indices;
}

template<typename indexType,size_t T>
indexType sub2ind(const std::array<indexType,T> & indices, const std::array<indexType,T> & dims)
{


    indexType idx=0;
    indexType multiplier=1;
    for(size_t i=0; i<T; ++i)
    {
        idx+=indices[i]*multiplier;
        multiplier*=dims[i];
    }
    return idx;

}

template<typename accessType, typename accessType2,typename dimType>
typename dimType::value_type sub2ind(const accessType & indices, const accessType2 &order, const dimType & dims)
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

    return idx;

}

/*! @} */
/*! @} */


}

#endif // INDEXUTIL_H

