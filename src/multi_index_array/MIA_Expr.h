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
#include <boost/mpl/at.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/comparison.hpp>

#include "kennytm/vtmp.hpp"
#include "Util.h"
#include "IndexUtil.h"
#include "ExprUtil.h"
#include "Index.h"
namespace LibMIA
{







namespace internal
{

/* beginning of taken from KennyTM's answer to this question:
http://stackoverflow.com/questions/10830406/copy-an-mplvector-c-to-a-static-array-at-compile-time*/
template <typename MPLVectorType>
class to_std_array
{
    typedef typename MPLVectorType::value_type element_type;
    static constexpr size_t length = boost::mpl::size<MPLVectorType>::value;
    typedef std::array<element_type, length> array_type;

    template <size_t... indices>
    static constexpr array_type
            make(const utils::vtmp::integers<indices...>&) noexcept
    {
        return array_type{{
            boost::mpl::at_c<MPLVectorType, indices>::type::value...
        }};
    }

public:
    static constexpr array_type make() noexcept
    {
        return make(utils::vtmp::iota<length>{});
    }
};
/*end of taken from KennyTM's answer to this question:
http://stackoverflow.com/questions/10830406/copy-an-mplvector-c-to-a-static-array-at-compile-time*/

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


template<typename _MIA, typename array1Type,typename array2Type, typename array3Type, bool hasOwnership>
struct lattice_maker
{

};

template<typename _MIA, typename array1Type,typename array2Type, typename array3Type>
struct lattice_maker<_MIA,array1Type,array2Type,array3Type,true>
{
    static auto apply(_MIA & _mia,const array1Type& array1, const array2Type& array2, const array3Type& array3)
    ->decltype(_mia.toLatticeDiscard(array1,array2,array3))
    {
        return _mia.toLatticeDiscard(array1,array2,array3);
    }
};

template<typename _MIA, typename array1Type,typename array2Type, typename array3Type>
struct lattice_maker<_MIA,array1Type,array2Type,array3Type,false>
{
    static auto apply(_MIA & _mia,const array1Type& array1, const array2Type& array2, const array3Type& array3)
    ->decltype(_mia.toLatticeExpression(array1,array2,array3))
    {
        return _mia.toLatticeExpression(array1,array2,array3);
    }
};

struct index_decrementer
{
    template<
        class _MIA,
        class Seq,
        bool hasownership,
        size_t inter_product_number,
        typename boost::enable_if<
            boost::mpl::greater<
                boost::mpl::size_t<inter_product_number>
                ,boost::mpl::size_t<0>
            >,
            int
        >::type=0
    >
    static auto apply(MIA_Atom<_MIA,Seq,hasownership,inter_product_number>& mia_atom)->MIA_Atom<_MIA,typename internal::decrement_back_indices<Seq,inter_product_number>::newSeq,hasownership>
    {
        _MIA * _mia=mia_atom.m_mia;
        mia_atom.m_mia=nullptr;
        return MIA_Atom<_MIA,typename internal::decrement_back_indices<Seq,inter_product_number>::newSeq,hasownership>(_mia);
    }

    template<
        class _MIA,
        class Seq,
        bool hasownership,
        size_t inter_product_number,
        typename boost::enable_if<
            boost::mpl::equal_to<
                boost::mpl::size_t<inter_product_number>
                ,boost::mpl::size_t<0>
            >,
        int >::type=0
    >
    static auto apply(MIA_Atom<_MIA,Seq,hasownership,inter_product_number>& mia_atom)->MIA_Atom<_MIA,Seq,hasownership,inter_product_number>
    {
        _MIA * _mia=mia_atom.m_mia;
        mia_atom.m_mia=nullptr;
        return MIA_Atom<_MIA,Seq,hasownership,inter_product_number>(_mia);
    }

};

template<class l_Seq,class r_Seq>
struct product_solve_expr_helper
{

    typedef typename internal::pull_left_operand_index_sequence<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value,0> left_index_order;
    typedef typename internal::pull_left_operand_index_sequence<r_Seq,l_Seq,boost::mpl::empty<r_Seq>::value,0> right_outer_index_order;
    typedef internal::pull_right_index_order<l_Seq,r_Seq,boost::mpl::empty<r_Seq>::value> pulling_index_order;


    template<class l_MIA_type, class r_MIA_type>
    static std::array<typename l_MIA_type::index_type,
                boost::mpl::size<typename left_index_order::no_match_order_sequence>::value +
                boost::mpl::size<typename right_outer_index_order::no_match_order_sequence>::value+
                boost::mpl::size<typename left_index_order::inter_match_order_sequence>::value>
        run(const l_MIA_type& l_MIA, const r_MIA_type& r_MIA)
    {
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<l_Seq,internal::binary_rule>::run();
        //check to makes sure no right-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();
        perform_cartesian_check<l_Seq,r_Seq,internal::product_rule>::run();

        //        print_array(left_outer_product_order, "left_outer");
        //        print_array(left_inner_product_order,"left_inner");
        //        print_array(left_inter_product_order,"left_inter");
        //
        //        print_array(right_outer_product_order, "right_outer");
        //        print_array(right_inner_product_order,"right_inner");
        //        print_array(right_inter_product_order,"right_inter");
        std::array<typename l_MIA_type::index_type,
                boost::mpl::size<typename left_index_order::no_match_order_sequence>::value +
                boost::mpl::size<typename right_outer_index_order::no_match_order_sequence>::value+
                boost::mpl::size<typename left_index_order::inter_match_order_sequence>::value> cMIA_dims;
        size_t curIdx=0;
        internal::reorder_from(l_MIA.dims(),left_outer_product_order(),cMIA_dims,curIdx);
        internal::reorder_from(r_MIA.dims(),right_outer_product_order(),cMIA_dims,curIdx);
        internal::reorder_from(l_MIA.dims(),internal::to_std_array<typename left_index_order::inter_match_order_sequence>::make(),cMIA_dims,curIdx);
        return cMIA_dims;

    }
    static constexpr auto left_outer_product_order()->decltype(internal::to_std_array<typename left_index_order::no_match_order_sequence>::make())
    {
        return internal::to_std_array<typename left_index_order::no_match_order_sequence>::make();
    }
    static constexpr auto right_outer_product_order()->decltype(internal::to_std_array<typename right_outer_index_order::no_match_order_sequence>::make())
    {
        return internal::to_std_array<typename right_outer_index_order::no_match_order_sequence>::make();
    }
    static constexpr auto left_inner_product_order()->decltype(internal::to_std_array<typename left_index_order::match_order_sequence>::make())
    {
        return internal::to_std_array<typename left_index_order::match_order_sequence>::make();
    }
    static constexpr auto right_inner_product_order()->decltype(internal::to_std_array<typename pulling_index_order::match_order>::make())
    {
        return internal::to_std_array<typename pulling_index_order::match_order>::make();
    }
    static constexpr auto left_inter_product_order()->decltype(internal::to_std_array<typename left_index_order::inter_match_order_sequence>::make())
    {
        return internal::to_std_array<typename left_index_order::inter_match_order_sequence>::make();
    }
    static constexpr auto right_inter_product_order()->decltype(internal::to_std_array<typename pulling_index_order::inter_match_order>::make())
    {
        return internal::to_std_array<typename pulling_index_order::inter_match_order>::make();
    }
};








template<class _MIA,class m_Seq,bool ownership,size_t inter_product_number>
class MIA_Atom
{
public:
    _MIA* m_mia;
    static constexpr bool mHasOwnership=ownership;
    MIA_Atom(_MIA* mia): m_mia(mia)
    {
        static_assert(internal::is_MIA<_MIA>::value,"Somehow expression was instantiated with a non-MIA class.");


    }
    ~MIA_Atom(){
        if(m_mia && mHasOwnership)
            delete m_mia;

    }

    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    auto operator*(const MIA_Atom<otherMIA,r_Seq, other_ownership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            true,
            MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    {



        return perform_product<otherMIA,r_Seq,other_ownership,other_inter_number>(Rhs);


    }

    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    auto operator|(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            true,
            MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    {

        typedef product_solve_expr_helper<m_Seq,r_Seq> helper;

        auto cMIA_dims=helper::run(*m_mia,*(Rhs.m_mia));


        auto aLat=m_mia->toLatticeExpression(helper::left_inner_product_order(),helper::left_outer_product_order(),helper::left_inter_product_order());



        auto bLat=Rhs.m_mia->toLatticeExpression(helper::right_inner_product_order(),helper::right_outer_product_order(),helper::right_inter_product_order());

        auto cLat=aLat.solve(bLat);



        typedef typename MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        constexpr size_t inter_product_number=MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
        MIA_return_type* cMIA(new MIA_return_type(cMIA_dims,std::move(cLat)));


        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,inter_product_number>(cMIA);




    }




    MIA_Atom& operator=(const MIA_Atom & Rhs)
    {


        if(Rhs.mHasOwnership)
            *m_mia=std::move(*(Rhs.m_mia));
        else
            *m_mia=*(Rhs.m_mia);
        return *this;
    }

    template<class otherMIA,bool other_ownership,size_t other_inter_number>
    MIA_Atom& operator=(const MIA_Atom<otherMIA,m_Seq,other_ownership,other_inter_number> & Rhs)
    {



        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform assignment.");

        if(other_ownership)
            *m_mia=std::move(*(Rhs.m_mia));
        else
            *m_mia=*(Rhs.m_mia);
        return *this;



    }

    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    MIA_Atom& operator=(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)
    {


        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform assignment.");


        typedef perform_cartesian_check<m_Seq,r_Seq,internal::assign_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();

        typedef internal::pull_match_order<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;


        if(other_ownership){

            m_mia->assign(std::move(*(Rhs.m_mia)),internal::to_std_array<typename pulling_index_order::match_order>::make());
        }
        else
            m_mia->assign(*(Rhs.m_mia),internal::to_std_array<typename pulling_index_order::match_order>::make());


        return *this;


    }

    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    MIA_Atom& operator+=(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)
    {


        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform addition.");

        typedef perform_cartesian_check<m_Seq,r_Seq,internal::merge_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();


        typedef internal::pull_match_order<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;


        m_mia->plus_equal(*(Rhs.m_mia),internal::to_std_array<typename pulling_index_order::match_order>::make());

        return *this;


    }


    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    auto operator+(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAMergeReturnType<_MIA,otherMIA>::type,
            m_Seq,
            true,
            inter_product_number
        >
    {



        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform addition.");

        typedef perform_cartesian_check<m_Seq,r_Seq,internal::merge_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();
        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;
        typedef internal::pull_match_order<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;



        cType* cMIA(new cType(m_mia->plus_(*(Rhs.m_mia),internal::to_std_array<typename pulling_index_order::match_order>::make())));

        MIA_Atom<cType,m_Seq,true,inter_product_number> C(cMIA);
        return C;


    }


    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    MIA_Atom& operator-=(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)
    {


        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform addition.");

        typedef perform_cartesian_check<m_Seq,r_Seq,internal::merge_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();

        typedef internal::pull_match_order<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;



        m_mia->minus_equal(*(Rhs.m_mia),internal::to_std_array<typename pulling_index_order::match_order>::make());

        return *this;


    }


    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    auto operator-(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAMergeReturnType<_MIA,otherMIA>::type,
            m_Seq,
            true,
            inter_product_number
        >
    {


        static_assert(internal::order<_MIA>::value==internal::order<otherMIA>::value,"Orders of two MIAs must be the same to perform addition.");

        typedef perform_cartesian_check<m_Seq,r_Seq,internal::merge_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<m_Seq,internal::binary_rule>::run();
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();
        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;
        typedef internal::pull_match_order<m_Seq,r_Seq,boost::mpl::empty<m_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;

        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;
        cType* cMIA=new cType(m_mia->minus_(*(Rhs.m_mia),internal::to_std_array<typename pulling_index_order::match_order>::make()));
        MIA_Atom<cType,m_Seq,true,inter_product_number> C(cMIA);
        return C;


    }


    auto operator~()->MIA_Atom<_MIA,typename internal::decrement_back_indices<m_Seq,inter_product_number>::newSeq,mHasOwnership>{
        return index_decrementer::template apply<_MIA,m_Seq,mHasOwnership,inter_product_number>(*this);

    }





private:




    template<class otherMIA,class r_Seq,bool otherOwnership,size_t other_inter_number,
                typename boost::disable_if_c<isPureOuterInterProduct<m_Seq,r_Seq>::type::value,int>::type=0
            >
    auto perform_product(const MIA_Atom<otherMIA,r_Seq, otherOwnership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            true,
            MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    {



        //std::cout << "Entered MIA Expr*" << std::endl;
        typedef product_solve_expr_helper<m_Seq,r_Seq> helper;

        auto cMIA_dims=helper::run(*m_mia,*(Rhs.m_mia));


        auto aLat=lattice_maker<_MIA,decltype(helper::left_outer_product_order()),decltype(helper::left_inner_product_order()),decltype(helper::left_inter_product_order()),mHasOwnership>
            ::apply(*m_mia,helper::left_outer_product_order(),helper::left_inner_product_order(),helper::left_inter_product_order());

        auto bLat=lattice_maker<otherMIA,decltype(helper::right_inner_product_order()),decltype(helper::right_outer_product_order()),decltype(helper::right_inter_product_order()),otherOwnership>
            ::apply(*(Rhs.m_mia),helper::right_inner_product_order(),helper::right_outer_product_order(),helper::right_inter_product_order());

       // std::cout << "ALat done" << std::endl;
        //auto bLat=Rhs.m_mia->toLatticeExpression(helper::right_inner_product_order(),helper::right_outer_product_order(),helper::right_inter_product_order());
        //std::cout << "BLat done" << std::endl;
        auto cLat=aLat*bLat;
        //std::cout << "CLat done" << std::endl;
        //cLat.print();



        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        constexpr size_t inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
        MIA_return_type* cMIA(new MIA_return_type(cMIA_dims,std::move(cLat)));
        //std::cout << "Lattice product finished " << std::endl;

        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,inter_product_number>(cMIA);




    }
    //when multiplication only has inter or outer products, we don't need to create a lattice, call specialized function in mia
    template<class otherMIA,class r_Seq,size_t other_inter_number,typename boost::enable_if_c<isPureOuterInterProduct<m_Seq,r_Seq>::type::value,int>::type=0>
    auto perform_product(const MIA_Atom<otherMIA,r_Seq, other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            true,
            MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    {


        //std::cout << "Pure outer/inter started " << std::endl;
        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        typedef product_solve_expr_helper<m_Seq,r_Seq> helper;


        MIA_return_type* cMIA(new MIA_return_type(m_mia->noLatticeMult(*Rhs.m_mia,helper::left_inter_product_order(),
                                                                       helper::left_outer_product_order(),
                                                                       helper::right_inter_product_order(),
                                                                       helper::right_outer_product_order())));
        constexpr size_t inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;

        //std::cout << "Pure outer/inter finished " << std::endl;
        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,inter_product_number>(cMIA);




    }

};





}// namespace LibMIA

#endif // MIA_EXPR_H_INCLUDED
