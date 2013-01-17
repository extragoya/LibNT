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

#define __UNIQUE_ID__ __COUNTER__           //should work in gcc 4.3 and later and also MS compiler
#define MIAINDEX ProdInd<__UNIQUE_ID__>
namespace LibMIA
{



template<char i>
struct ElemWiseInd {};



template<size_t id,bool ElemWise=false>
struct ProdInd
{
    constexpr static bool elemval=ElemWise;
    ProdInd<id,true> operator !()
    {
        return ProdInd<id,!ElemWise>();
    }

};



//these provide grammar rules for MIA expressions
namespace internal
{

constexpr int product_rule=0;
constexpr int assign_rule=1;
constexpr int binary_rule=0;

template<class T>
struct is_ProdInd: public boost::false_type {};

template<size_t i>
struct is_ProdInd<ProdInd<i>>: public boost::true_type {};

template<class T1,int T2>
struct match_rule;

//number of matches allowed with an ordinary index during an MIA product
template<size_t label>
struct match_rule<ProdInd<label,false>,product_rule>
{
    typedef boost::mpl::vector_c<int,0,1> allowed_matches;

};

//number of matches allowed with an elemwise index during an MIA product
template<size_t label>
struct match_rule<ProdInd<label,true>,product_rule>
{
    typedef boost::mpl::vector_c<int,1> allowed_matches;

};

//number of matches allowed with an index during an MIA assignment
template<size_t label,bool elemwise>
struct match_rule<ProdInd<label,elemwise>,assign_rule>
{
    typedef boost::mpl::vector_c<int,1> allowed_matches;

};



template<class T1,int T2>
struct auto_match_rule;

//use when we want to ensure no MIA has a repeated index within the same MIA
template<size_t id,bool elemwise>
struct auto_match_rule<ProdInd<id,elemwise>,binary_rule>
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


} //namespace internal




} //namespace LibMIA


#endif // INDEX_H_INCLUDED
