// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef INDEX_H_INCLUDED
#define INDEX_H_INCLUDED

#include <string>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/min_max.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/integral_c_tag.hpp>
#include <boost/mpl/aux_/config/static_constant.hpp>
#define __UNIQUE_ID__ __COUNTER__           //should work in gcc 4.3 and later and also MS compiler
#define MIAINDEX ProdInd<__UNIQUE_ID__>
namespace LibMIA
{







template<size_t ID,int ElemWise=0>
struct ProdInd
{
    constexpr static int elemval=ElemWise;
    constexpr static size_t id=ID;
    typedef typename
        boost::mpl::max<
            boost::mpl::int_<elemval-1>,
            boost::mpl::int_<0>
        >::type max_elem;
    typedef ProdInd<
        id,
        max_elem::value
    > decrement_type;
    ProdInd<id,ElemWise+1> operator !()
    {
        return ProdInd<ID,ElemWise+1>();
    }

};



//these provide grammar rules for MIA expressions
namespace internal
{

constexpr int product_rule=0;
constexpr int assign_rule=1;
constexpr int merge_rule=2;
constexpr int binary_rule=0;

template<class T>
struct is_ProdInd: public boost::false_type {};

template<size_t i,int elemval>
struct is_ProdInd<ProdInd<i,elemval>>: public boost::true_type {};

template<class T,int rule_id>
struct match_rule;

//number of matches allowed with an elemwise index during an MIA product
template<size_t id,int elemval>
struct match_rule<ProdInd<id,elemval>,product_rule>
{
    typedef boost::mpl::vector_c<int,1> allowed_matches;

};

//number of matches allowed with an ordinary index during an MIA product
template<size_t id>
struct match_rule<ProdInd<id,0>,product_rule>
{
    typedef boost::mpl::vector_c<int,0,1> allowed_matches;

};



//number of matches allowed with an index during an MIA assignment
template<size_t id,int elemval>
struct match_rule<ProdInd<id,elemval>,assign_rule>
{
    typedef boost::mpl::vector_c<int,1> allowed_matches;

};


//number of matches allowed with an ordinary index during an MIA merge (add/subtract)
template<size_t id,int elemval>
struct match_rule<ProdInd<id,elemval>,merge_rule>
{
    typedef boost::mpl::vector_c<int,1> allowed_matches;

};


template<class T,int rule_id>
struct auto_match_rule;

//use when we want to ensure no MIA has a repeated index within the same MIA
template<size_t id,int elemval>
struct auto_match_rule<ProdInd<id,elemval>,binary_rule>
{
    typedef boost::mpl::vector_c<int,0> allowed_matches;

};

//******BEGIN taken from
//http://blog.shandyba.com/2009/12/17/converting-variadic-template-arguments-pack-to-boost-mpl-sequence/
//March 26 2012
//*****************************************
//General definition of the helper class
template <typename ...Args> struct FromVariadic;

//This specialization does the actual job: it splits the whole pack
//into 2 parts: one single type T and the rest of types Args...
//As soon as it is done T is added to an mpl::vector.
//"bottom--up" recursion is used to fetch all types
template <class T, typename ...Args>
struct FromVariadic<T, Args...>
{
    typedef typename boost::mpl::push_front<typename FromVariadic<Args...>::type, T>::type type;
};

//This is a specialization for the case when only one type is passed
//and also finishes recursive descent
template <class T>
struct FromVariadic<T>
{
    typedef boost::mpl::vector<T> type;
};

//This one handles the case when no types where passed at all
template <>
struct FromVariadic<>
{
    typedef boost::mpl::vector<> type;
};
//******END taken from
//http://blog.shandyba.com/2009/12/17/converting-variadic-template-arguments-pack-to-boost-mpl-sequence/
//March 26 2012
//*****************************************

template<class ...Ts>
struct Indicial_Sequence
{
    typedef typename internal::FromVariadic<Ts...>::type sequence;
};

template<typename T1, typename T2>
struct same_product_index_id: public boost::false_type {};

template<size_t id1, int elem_val1,size_t id2, int elem_val2>
struct same_product_index_id<ProdInd<id1,elem_val1>,ProdInd<id2,elem_val2>> {

    BOOST_STATIC_CONSTANT(bool, value = id1==id2);
    typedef boost::mpl::integral_c_tag tag;
    typedef same_product_index_id type;
    typedef bool value_type;
    operator bool() const { return this->value; }
};



} //namespace internal




} //namespace LibMIA


#endif // INDEX_H_INCLUDED
