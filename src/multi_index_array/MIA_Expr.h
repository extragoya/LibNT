#ifndef MIA_EXPR_H_INCLUDED
#define MIA_EXPR_H_INCLUDED

#include<array>
#include <type_traits>
#include <tuple>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/for_each.hpp>


#include "Util.h"
#include "Index.h"
namespace LibMIA
{



namespace internal
{


template< typename Sequence >
class sequence_array : public std::array< typename Sequence::value_type, boost::mpl::size<Sequence>::type::value>
{
    typedef typename std::array< typename Sequence::value_type, boost::mpl::size<Sequence>::type::value>::iterator iterator;
    struct copier_
    {
        copier_(iterator it) : it_(it) {}
        template<typename U> void operator()(U u)
        {
            *(it_++) = u;
        }
        iterator it_;
    };
public:
    sequence_array()
    {
        boost::mpl::for_each<Sequence>(copier_(this->begin()));
    }
};

//operator*(MIA_Expression)

};

template<class _MIA,class m_Seq>
class MIA_Atom
{
public:

    MIA_Atom(_MIA& mia): m_mia(mia)
    {
        static_assert(internal::is_MIA<_MIA>::value,"Somehow expression was instantiated with a non-MIA class.");


    }

    template<class otherMIA,class r_Seq>
    auto operator*(const MIA_Atom<otherMIA,r_Seq> & Rhs)->MIA_Atom<typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence>
    {

        //typedef typename internal::Indicial_Sequence<Rhs_inds...>::sequence r_Seq;


        typedef typename internal::cartesian_product_indices<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,internal::product_rule,0> left_sequence_check;
        typedef typename internal::cartesian_product_indices<r_Seq,m_Seq,boost::mpl::empty<r_Seq>::value,internal::product_rule,0> right_sequence_check;

        //check to makes sure no left-hand indice is repeated
        static_assert(internal::auto_cartesian_product_indices<m_Seq,boost::mpl::empty<m_Seq>::value,internal::binary_rule>::value,"Repeated index in left-hand operand.");
        //check to make sure no right-hand indice is repeated
        static_assert(internal::auto_cartesian_product_indices<r_Seq,boost::mpl::empty<r_Seq>::value,internal::binary_rule>::value,"Repeated index in right-hand operand.");
        //now check to make sure left and right-hand operands match up properly
        static_assert(left_sequence_check::value,"A left-hand index does not match up properly with a right-hand index.");
        static_assert(right_sequence_check::value,"A right-hand index does not match up properly with a left-hand index.");

        internal::sequence_array<typename left_sequence_check::match_order_sequence> left_inner_product_order;
        internal::sequence_array<typename left_sequence_check::no_match_order_sequence> left_outer_product_order;
        internal::sequence_array<typename left_sequence_check::inter_match_order_sequence> left_inter_product_order;

        internal::sequence_array<typename right_sequence_check::match_order_sequence> right_inner_product_order;
        internal::sequence_array<typename right_sequence_check::no_match_order_sequence> right_outer_product_order;
        internal::sequence_array<typename right_sequence_check::inter_match_order_sequence> right_inter_product_order;

        //perform_mult(Rhs.m_mia,left_inner_product_order,left_outer_product_order,left_inter_product_order,right_inner_product_order,right_outer_product_order,right_inter_product_order);
        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        MIA_return_type test(5,5);

        auto cLat=m_mia.toLatticeExpression(left_outer_product_order,left_inner_product_order,left_inter_product_order)*Rhs.m_mia.toLatticeExpression(left_inner_product_order,left_outer_product_order,left_inter_product_order);
        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence>(test);



    }

    template<class otherMIA,class r_Seq>
    MIA_Atom & operator=(const MIA_Atom<otherMIA,r_Seq> & Rhs){


    }




    _MIA& m_mia;

private:




    //convert variadic templates to a boost::mpl sequence so that they can be worked with
    //typedef typename internal::FromVariadic<Ts...>::type m_Seq;
    //std::tuple<Ts...> m_tuple;

};

}// namespace LibMIA

#endif // MIA_EXPR_H_INCLUDED
