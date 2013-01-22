// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


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
#include "IndexUtil.h"
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

}

template<class l_Seq,class r_Seq,int rule>
struct perform_cartesian_check
{
    typedef typename internal::cartesian_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value,rule,0> left_sequence_check;
    typedef typename internal::cartesian_product_indices<r_Seq,l_Seq,boost::mpl::empty<r_Seq>::value,rule,0> right_sequence_check;
    static void run()
    {
        static_assert(left_sequence_check::value,"A left-hand index does not match up properly with a right-hand index.");
        static_assert(right_sequence_check::value,"A right-hand index does not match up properly with a left-hand index.");
    }

};

template<class Seq, int rule>
struct perform_auto_check
{
    typedef internal::auto_cartesian_product_indices<Seq,boost::mpl::empty<Seq>::value,rule> sequence_check;
    static void run()
    {
        static_assert(sequence_check::value,"Repeated index in operand.");
    }

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



        auto cLat=m_mia.toLatticeExpression(left_outer_product_order,left_inner_product_order,left_inter_product_order)*Rhs.m_mia.toLatticeExpression(left_inner_product_order,left_outer_product_order,left_inter_product_order);
        std::array<typename _MIA::index_type,left_outer_product_order.size()+right_outer_product_order.size()+left_inter_product_order.size()> cMIA_dims;
        size_t curIdx=0;
        internal::collect_dimensions(m_mia,left_outer_product_order,cMIA_dims,curIdx);
        internal::collect_dimensions(Rhs.m_mia,right_outer_product_order,cMIA_dims,curIdx);
        internal::collect_dimensions(m_mia,left_inter_product_order,cMIA_dims,curIdx);
        MIA_return_type cMIA(cMIA_dims,cLat.release_memptr());
        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence>(cMIA);



    }


    MIA_Atom& operator=(const MIA_Atom & Rhs)
    {


        //TODO check that no index is an inter product and check that orders match


        m_mia=Rhs.m_mia;



    }

    template<class otherMIA>
    MIA_Atom& operator=(const MIA_Atom<otherMIA,m_Seq> & Rhs)
    {

        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform assignment.");
        //TODO check that no index is an inter product and check that orders match


        m_mia=Rhs.m_mia;



    }

    template<class otherMIA,class r_Seq>
    MIA_Atom& operator=(const MIA_Atom<otherMIA,r_Seq> & Rhs)
    {

        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform assignment.");
        //TODO check that no index is an inter product and check that orders match
        typedef perform_cartesian_check<m_Seq,r_Seq,internal::assign_rule> cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();

        typedef internal::pull_right_index_order<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value> pulling_index_order;

        internal::sequence_array<typename pulling_index_order::match_order> right_assignment_order;
        std::array<typename otherMIA::index_type,internal::order<otherMIA>::value> r_Dims;

        m_mia.assign(Rhs.m_mia,right_assignment_order);



    }









    _MIA& m_mia;

private:





    //convert variadic templates to a boost::mpl sequence so that they can be worked with
    //typedef typename internal::FromVariadic<Ts...>::type m_Seq;
    //std::tuple<Ts...> m_tuple;

};

}// namespace LibMIA

#endif // MIA_EXPR_H_INCLUDED
