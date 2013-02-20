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
#include <boost/mpl/quote.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/comparison.hpp>

#include "Util.h"
#include "IndexUtil.h"
#include "ExprUtil.h"
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

template<class l_Seq,class r_Seq,int rule,typename Pred=boost::mpl::quote2<boost::is_same> >
struct perform_cartesian_check
{
    typedef typename internal::check_cartesian_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value,rule,Pred> left_sequence_check;
    typedef typename internal::check_cartesian_product_indices<r_Seq,l_Seq,boost::mpl::empty<r_Seq>::value,rule,Pred> right_sequence_check;
    static void run()
    {
        static_assert(left_sequence_check::value,"A left-hand index does not match up properly with a right-hand index.");
        static_assert(right_sequence_check::value,"A right-hand index does not match up properly with a left-hand index.");
    }

};

template<class Seq, int rule,typename Pred=boost::mpl::quote2<internal::same_product_index_id> >
struct perform_auto_check
{
    typedef internal::auto_cartesian_product_indices<Seq,boost::mpl::empty<Seq>::value,rule,Pred> sequence_check;
    static void run()
    {
        static_assert(sequence_check::value,"Repeated index in operand.");
    }

};


struct index_decrementer
{
    template<
        class _MIA,
        class Seq,
        size_t inter_product_number,
        typename boost::enable_if<
            boost::mpl::greater<
                boost::mpl::size_t<inter_product_number>
                ,boost::mpl::size_t<0>
            >,
            int
        >::type=0
    >
    static auto apply(MIA_Atom<_MIA,Seq,inter_product_number>& mia_atom)->MIA_Atom<_MIA,typename internal::decrement_back_indices<Seq,inter_product_number>::newSeq>
    {
        bool old_has_ownership=mia_atom.mHasOwnership;
        mia_atom.mHasOwnership=false;
        return MIA_Atom<_MIA,typename internal::decrement_back_indices<Seq,inter_product_number>::newSeq>(mia_atom.m_mia,old_has_ownership);
    }

    template<
        class _MIA,
        class Seq,
        size_t inter_product_number,
        typename boost::enable_if<
            boost::mpl::equal_to<
                boost::mpl::size_t<inter_product_number>
                ,boost::mpl::size_t<0>
            >,
        int >::type=0
    >
    static auto apply(MIA_Atom<_MIA,Seq,inter_product_number>& mia_atom)->MIA_Atom<_MIA,typename internal::decrement_back_indices<Seq,inter_product_number>::newSeq>
    {
        return mia_atom;
    }

};



template<class _MIA,class m_Seq,size_t inter_product_number>
class MIA_Atom
{
public:

    bool mHasOwnership;
    MIA_Atom(_MIA* mia,bool _has_ownership=false): m_mia(mia),mHasOwnership(_has_ownership)
    {
        static_assert(internal::is_MIA<_MIA>::value,"Somehow expression was instantiated with a non-MIA class.");


    }
    ~MIA_Atom(){
        if(mHasOwnership)
            delete m_mia;

    }

    template<class otherMIA,class r_Seq,size_t other_inter_number>
    auto operator*(const MIA_Atom<otherMIA,r_Seq, other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    { return 0;

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();
        perform_cartesian_check<m_Seq,r_Seq,internal::product_rule>::run();

        typedef typename internal::pull_left_operand_index_sequence<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,0> left_index_order;
        typedef typename internal::pull_left_operand_index_sequence<r_Seq,m_Seq,boost::mpl::empty<r_Seq>::value,0> right_outer_index_order;
        typedef internal::pull_right_index_order<m_Seq,r_Seq,boost::mpl::empty<r_Seq>::value> pulling_index_order;

        internal::sequence_array<typename left_index_order::match_order_sequence> left_inner_product_order;
        internal::sequence_array<typename left_index_order::no_match_order_sequence> left_outer_product_order;
        internal::sequence_array<typename left_index_order::inter_match_order_sequence> left_inter_product_order;

        internal::sequence_array<typename pulling_index_order::match_order> right_inner_product_order;
        internal::sequence_array<typename right_outer_index_order::no_match_order_sequence> right_outer_product_order;
        internal::sequence_array<typename pulling_index_order::inter_match_order> right_inter_product_order;




//        print_array(left_outer_product_order, "left_outer");
//        print_array(left_inner_product_order,"left_inner");
//        print_array(left_inter_product_order,"left_inter");
//
//        print_array(right_outer_product_order, "right_outer");
//        print_array(right_inner_product_order,"right_inner");
//        print_array(right_inter_product_order,"right_inter");

        auto cLat=m_mia->toLatticeExpression(left_outer_product_order,left_inner_product_order,left_inter_product_order)*Rhs.m_mia->toLatticeExpression(right_inner_product_order,right_outer_product_order,right_inter_product_order);
        //std::cout << "cLat " << std::endl;
        //cLat.print();
        std::array<typename _MIA::index_type,left_outer_product_order.size()+right_outer_product_order.size()+left_inter_product_order.size()> cMIA_dims;
        size_t curIdx=0;
        internal::collect_dimensions_from_order(*m_mia,left_outer_product_order,cMIA_dims,curIdx);
        internal::collect_dimensions_from_order(*(Rhs.m_mia),right_outer_product_order,cMIA_dims,curIdx);
        internal::collect_dimensions_from_order(*m_mia,left_inter_product_order,cMIA_dims,curIdx);

        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        constexpr size_t inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
        MIA_return_type* cMIA(new MIA_return_type(cMIA_dims,cLat.release_memptr()));


        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,inter_product_number>(cMIA,true);




    }

//    template<class otherMIA,class r_Seq>
//    auto operator|(const MIA_Atom<otherMIA,r_Seq> & Rhs)->MIA_Atom<typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number>
//    {
//
//        //check to makes sure no left-hand indice is repeated
//        perform_auto_check<m_Seq,internal::binary_rule>::run();
//        //check to makes sure no left-hand indice is repeated
//        perform_auto_check<r_Seq,internal::binary_rule>::run();
//        perform_cartesian_check<m_Seq,r_Seq,internal::product_rule>::run();
//
//        typedef typename internal::pull_left_operand_index_sequence<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,0> left_index_order;
//        typedef typename internal::pull_left_operand_index_sequence<r_Seq,m_Seq,boost::mpl::empty<r_Seq>::value,0> right_outer_index_order;
//        typedef internal::pull_right_index_order<m_Seq,r_Seq,boost::mpl::empty<r_Seq>::value> pulling_index_order;
//
//        internal::sequence_array<typename left_index_order::match_order_sequence> left_inner_product_order;
//        internal::sequence_array<typename left_index_order::no_match_order_sequence> left_outer_product_order;
//        internal::sequence_array<typename left_index_order::inter_match_order_sequence> left_inter_product_order;
//
//        internal::sequence_array<typename pulling_index_order::match_order> right_inner_product_order;
//        internal::sequence_array<typename right_outer_index_order::no_match_order_sequence> right_outer_product_order;
//        internal::sequence_array<typename pulling_index_order::inter_match_order> right_inter_product_order;
//
//
//
//
////        print_array(left_outer_product_order, "left_outer");
////        print_array(left_inner_product_order,"left_inner");
////        print_array(left_inter_product_order,"left_inter");
////
////        print_array(right_outer_product_order, "right_outer");
////        print_array(right_inner_product_order,"right_inner");
////        print_array(right_inter_product_order,"right_inter");
//
//        auto cLat=m_mia->toLatticeExpression(left_outer_product_order,left_inner_product_order,left_inter_product_order)*Rhs.m_mia->toLatticeExpression(right_inner_product_order,right_outer_product_order,right_inter_product_order);
//        //std::cout << "cLat " << std::endl;
//        //cLat.print();
//        std::array<typename _MIA::index_type,left_outer_product_order.size()+right_outer_product_order.size()+left_inter_product_order.size()> cMIA_dims;
//        size_t curIdx=0;
//        internal::collect_dimensions_from_order(*m_mia,left_outer_product_order,cMIA_dims,curIdx);
//        internal::collect_dimensions_from_order(*(Rhs.m_mia),right_outer_product_order,cMIA_dims,curIdx);
//        internal::collect_dimensions_from_order(*m_mia,left_inter_product_order,cMIA_dims,curIdx);
//
//        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
//        constexpr size_t inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
//        MIA_return_type* cMIA(new MIA_return_type(cMIA_dims,cLat.release_memptr()));
//
//
//        //create an MIA from cLat
//        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,inter_product_number>(cMIA,true);
//
//
//
//
//    }




    MIA_Atom& operator=(const MIA_Atom & Rhs)
    {


        *m_mia=*(Rhs.m_mia);
    }

    template<class otherMIA,size_t other_inter_number>
    MIA_Atom& operator=(const MIA_Atom<otherMIA,m_Seq,other_inter_number> & Rhs)
    {


        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform assignment.");
        *(m_mia)=*(Rhs.m_mia);



    }

    template<class otherMIA,class r_Seq,size_t other_inter_number>
    MIA_Atom& operator=(const MIA_Atom<otherMIA,r_Seq,other_inter_number> & Rhs)
    {


        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform assignment.");
        //TODO check that no index is an inter product and check that orders match
        typedef perform_cartesian_check<m_Seq,r_Seq,internal::assign_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();

        typedef internal::pull_match_order<r_Seq,m_Seq,boost::mpl::empty<r_Seq>::value,boost::mpl::quote2<internal::same_product_index_id> > pulling_index_order;

        internal::sequence_array<typename pulling_index_order::match_order> right_assignment_order;
        //print_array(right_assignment_order, "assign_order");
        std::array<typename otherMIA::index_type,internal::order<otherMIA>::value> r_Dims;

        m_mia->assign(*(Rhs.m_mia),right_assignment_order);

        return *this;


    }





    auto operator~()->MIA_Atom<_MIA,typename internal::decrement_back_indices<m_Seq,inter_product_number>::newSeq>{
        return index_decrementer::template apply<_MIA,m_Seq,inter_product_number>(*this);

    }



    _MIA* m_mia;

private:



};





}// namespace LibMIA

#endif // MIA_EXPR_H_INCLUDED
