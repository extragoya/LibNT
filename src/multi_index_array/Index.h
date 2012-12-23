#ifndef INDEX_H_INCLUDED
#define INDEX_H_INCLUDED

#include <string>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector_c.hpp>

#define __UNIQUE_ID__ __COUNTER__           //should work in gcc 4.3 and later and also MS compiler
#define PRODINDEX ProdInd<__UNIQUE_ID__>
namespace LibMIA
{



template<char i>
struct ElemWiseInd {};



template<int id,bool ElemWise=false>
struct ProdInd
{
    constexpr static bool elemval=ElemWise;
    ProdInd<id,true> operator !()
    {
        return ProdInd<id,!ElemWise>();
    }

};




namespace internal
{

constexpr int product_rule=0;
constexpr int binary_rule=0;

template<class T>
struct is_ProdInd: public boost::false_type {};

template<int i>
struct is_ProdInd<ProdInd<i>>: public boost::true_type {};

template<int i>
struct is_ProdInd<ElemWiseInd<i>>: public boost::true_type {};

template<class T1,int T2>
struct match_rule;

template<int label>
struct match_rule<ProdInd<label,false>,product_rule>
{
    typedef boost::mpl::vector_c<int,0,1> allowed_matches;

};

template<int label>
struct match_rule<ProdInd<label,true>,product_rule>
{
    typedef boost::mpl::vector_c<int,1> allowed_matches;

};

template<class T1,int T2>
struct auto_match_rule;

template<int id,bool elemwise>
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
