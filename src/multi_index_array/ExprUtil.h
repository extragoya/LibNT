// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef EXPRUTIL_H
#define EXPRUTIL_H




#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/print.hpp>
#include <boost/mpl/count.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/print.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/bind.hpp>
#include <boost/mpl/erase.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/utility/enable_if.hpp>

#include "Index.h"
#include "Util.h"
using namespace boost::mpl::placeholders;
namespace LibMIA
{

/** \addtogroup util Utilities
 *  @{
*/

/** \addtogroup expression_util MIA Expression Utilities
 *  @{
 */






namespace internal
{

//used to get the order of inner, inter, and outer product indices for left operands
template<typename LSeq,typename RSeq,bool empty_LSeq,int recursive_depth,typename Pred=boost::mpl::quote2<boost::is_same> >
struct pull_left_operand_index_sequence
{

    typedef boost::mpl::vector_c<int> match_order_sequence;
    typedef boost::mpl::vector_c<int> inter_match_order_sequence;
    typedef boost::mpl::vector_c<int> no_match_order_sequence;
};

template<typename LSeq,typename RSeq,int recursive_depth,typename Pred >
struct pull_left_operand_index_sequence<LSeq,RSeq,false,recursive_depth,Pred>
{

    //Pull first index type off of left sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<LSeq>::type>::type first_LSeq;
    //count its occurences in the right sequence
    typedef typename boost::mpl::count_if<
        RSeq,
        typename boost::mpl::apply_wrap2<
            Pred,
            boost::mpl::_1,
            first_LSeq
        >
    >::type n;

    //pop the first type off of LSeq
    typedef typename boost::mpl::pop_front<LSeq>::type poppedLSeq;

    //obtain index sequence
    typedef  pull_left_operand_index_sequence<poppedLSeq,RSeq,boost::mpl::empty<poppedLSeq>::value,recursive_depth+1,Pred> next_pull_left_operand_index_sequence;


    typedef typename boost::mpl::if_<
        boost::mpl::and_<
            boost::mpl::equal_to<
                n,
                boost::mpl::int_<1>
            >,
            boost::mpl::equal_to<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
        >,
        typename boost::mpl::push_front<
            typename next_pull_left_operand_index_sequence::match_order_sequence,
            boost::mpl::int_<recursive_depth>
        >::type,
        typename next_pull_left_operand_index_sequence::match_order_sequence
    >::type match_order_sequence;


    typedef typename boost::mpl::if_<
        boost::mpl::and_<
            boost::mpl::equal_to<
                n,boost::mpl::int_<0>
            >,
            boost::mpl::equal_to<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
        >,
        typename boost::mpl::push_front<
            typename next_pull_left_operand_index_sequence::no_match_order_sequence,
            boost::mpl::int_<recursive_depth>>::type,
        typename next_pull_left_operand_index_sequence::no_match_order_sequence
    >::type no_match_order_sequence;

    typedef typename boost::mpl::if_<
        boost::mpl::and_<
            boost::mpl::equal_to<
                n,boost::mpl::int_<1>
            >,
            boost::mpl::greater<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
        >,
        typename boost::mpl::push_front<
            typename next_pull_left_operand_index_sequence::inter_match_order_sequence,
            boost::mpl::int_<recursive_depth>
        >::type,
        typename next_pull_left_operand_index_sequence::inter_match_order_sequence
    >::type inter_match_order_sequence;
};

//base case - set to true. This only occurs when LSeq is an empty mpl sequence
template<typename LSeq,typename RSeq,bool empty_LSeq,int expression_type, typename Pred=boost::mpl::quote2<boost::is_same> >
struct check_cartesian_product_indices
{
    typedef boost::true_type allowed_recursive;



};

//this structure checks how many times the first element of LSeq occurs in RSeq
//Once this is done, it then checks to make sure the frequency follows the expression_type rules
//finally it recurses to check the next element of LSeq
template<typename LSeq,typename RSeq,int expression_type,typename Pred>
struct check_cartesian_product_indices<LSeq,RSeq,false,expression_type,Pred>
{

    //Pull first index type off of left sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<LSeq>::type>::type first_LSeq;
    //count its occurences in the right sequence
    typedef typename boost::mpl::count_if<
        RSeq,
        typename boost::mpl::apply_wrap2<
            Pred,
            boost::mpl::_1,
            first_LSeq
        >
    >::type n;
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
    typedef  check_cartesian_product_indices<poppedLSeq,RSeq,boost::mpl::empty<poppedLSeq>::value,expression_type,Pred> next_check_cartesian_product_indices;


    //now check the remaining LSeq types using recursive calls
    typedef typename boost::mpl::and_<
    allowed,
    typename next_check_cartesian_product_indices::allowed_recursive
    > allowed_recursive;
    constexpr static bool value= allowed_recursive::value;



};




//base case - set to true. This only occurs when LSeq is an empty mpl sequence
template<typename Seq,bool empty_LSeq,int expression_type, typename Pred=boost::mpl::quote2<internal::same_product_index_id> >
struct auto_cartesian_product_indices
{
    typedef boost::true_type allowed_recursive;
};

//this structure checks how many times the first element of Seq occurs in the rest of the sequence.
//Once this is done, it then checks to make sure the frequency follows the expression_type rules
//finally it recurses to check the next element of Seq
template<typename Seq,int expression_type, typename Pred>
struct auto_cartesian_product_indices<Seq,false,expression_type,Pred>
{

    //Pull first index type off of the sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<Seq>::type>::type first_Seq;
    //pop the first type off of Seq
    typedef typename boost::mpl::pop_front<Seq>::type poppedSeq;
    //count its occurences in the sequence
    typedef typename boost::mpl::count_if<
        poppedSeq,
        typename boost::mpl::apply_wrap2<
            Pred,
            boost::mpl::_1,
            first_Seq
        >
    >::type n;
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

template<typename Seq, typename Seq_toMatch,bool empty_LSeq,typename Pred=boost::mpl::quote2<boost::is_same> >
struct pull_match_order
{
    typedef boost::mpl::vector_c<int> match_order;

};

template<typename Seq, typename Seq_toMatch,typename Pred>
struct pull_match_order<Seq,Seq_toMatch,false,Pred>
{

    //Pull first index type off of left sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<Seq>::type>::type first_Seq;
    //find its position in RSeq
    typedef typename boost::mpl::find_if<
        Seq_toMatch,
        typename boost::mpl::apply_wrap2<
            Pred,
            boost::mpl::_1,
            first_Seq
        >
    >::type find_pos;



    typedef typename boost::mpl::not_<
                        boost::is_same<
                            find_pos,
                            typename boost::mpl::end<Seq_toMatch>::type
                        >
                    >::type was_found;

    //pop the first type off of Seq
    typedef typename boost::mpl::pop_front<Seq>::type poppedSeq;
    //obtain next match for rest of Seq (must include Pred!!)
    typedef  pull_match_order<poppedSeq,Seq_toMatch,boost::mpl::empty<poppedSeq>::value,Pred> next_match_order;

    //if first_LSeq was found in RSeq and it's elementwise, add the match location to the runing inter_match_order
    typedef typename boost::mpl::if_<
        was_found,
        typename boost::mpl::push_front<
            typename next_match_order::match_order,
            typename boost::mpl::distance<
                typename boost::mpl::begin<Seq_toMatch>::type,
                find_pos
            >::type
        >::type,
        typename next_match_order::match_order
    >::type match_order;


};

template<typename LSeq, typename RSeq,bool empty_LSeq,typename Pred=boost::mpl::quote2<boost::is_same> >
struct pull_right_index_order
{

    typedef boost::mpl::vector_c<int> match_order;
    typedef boost::mpl::vector_c<int> inter_match_order;


};


template<typename LSeq, typename RSeq,typename Pred>
struct pull_right_index_order<LSeq,RSeq,false,Pred>
{

    //Pull first index type off of left sequence
    typedef typename boost::mpl::deref<typename boost::mpl::begin<LSeq>::type>::type first_LSeq;
    //find its position in RSeq
    typedef typename boost::mpl::find_if<
        RSeq,
        typename boost::mpl::apply_wrap2<
            Pred,
            boost::mpl::_1,
            first_LSeq
        >
    >::type find_pos;



    typedef typename boost::mpl::not_<
                        boost::is_same<
                            find_pos,
                            typename boost::mpl::end<RSeq>::type
                        >
                    >::type was_found;

    //pop the first type off of LSeq
    typedef typename boost::mpl::pop_front<LSeq>::type poppedLSeq;
    //obtain next match for rest of LSeq
    typedef  pull_right_index_order<poppedLSeq,RSeq,boost::mpl::empty<poppedLSeq>::value,Pred> next_pull_right_index_order;

    //if first_LSeq was found in RSeq and it's elementwise, add the match location to the runing inter_match_order
    typedef typename boost::mpl::if_<
        boost::mpl::and_<
            was_found,
            boost::mpl::greater<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
        >,
        typename boost::mpl::push_front<
            typename next_pull_right_index_order::inter_match_order,
            typename boost::mpl::distance<
                typename boost::mpl::begin<RSeq>::type,
                find_pos
            >::type
        >::type,
        typename next_pull_right_index_order::inter_match_order
    >::type inter_match_order;

    //if first_LSeq was found in RSeq and it's not elementwise, add the match location to the runing match_order
    typedef typename boost::mpl::if_<
        boost::mpl::and_<
            was_found,
            boost::mpl::equal_to<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
        >,
        typename boost::mpl::push_front<
            typename next_pull_right_index_order::match_order,
            typename boost::mpl::distance<
                typename boost::mpl::begin<RSeq>::type,
                find_pos
            >::type
        >::type,
        typename next_pull_right_index_order::match_order
    >::type match_order;


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
            boost::mpl::greater<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
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
            boost::mpl::equal_to<
                boost::mpl::int_<first_LSeq::elemval>,
                boost::mpl::int_<0>
            >
        >,
        typename boost::mpl::push_front<
            typename next_pull_product_indices::outer_product_indices,
            first_LSeq
        >::type,
        typename next_pull_product_indices::outer_product_indices
    >::type outer_product_indices;

};

//base case, just an empty vector
template<class Seq, size_t range_to_decrement,class Enable = void>
struct get_decremented_prod_index
{
    typedef boost::mpl::vector<> decremented_Seq;
};

//only enabled when range_to_decrement>0
template<class Seq, size_t range_to_decrement>
struct get_decremented_prod_index<
        Seq,
        range_to_decrement,
        typename boost::enable_if<
            boost::mpl::greater<
                boost::mpl::size_t<range_to_decrement>,
                boost::mpl::size_t<0>
            >
        >::type
    >
{
    //Pull first index type off of sequence
    typedef typename boost::mpl::deref<
        typename boost::mpl::prior<
            typename boost::mpl::end<Seq>::type
        >::type
    >::type first_Seq;
    //pop the last seq from the seq
    typedef typename boost::mpl::pop_back<Seq>::type popped_Seq;


    typedef typename boost::mpl::push_back<
        typename get_decremented_prod_index<popped_Seq,range_to_decrement-1>::decremented_Seq,
        typename first_Seq::decrement_type
    >::type decremented_Seq;
};

template<class Seq,size_t range_to_decrement>
struct decrement_back_indices
{



    typedef typename boost::mpl::advance<
        typename boost::mpl::end<Seq>::type,
        boost::mpl::int_<-1*static_cast<int>(range_to_decrement)>
    >::type erase_iterator;

    typedef typename boost::mpl::erase<
        Seq,
        erase_iterator,
        typename boost::mpl::end<Seq>::type
    >::type erased_Seq;
    typedef typename get_decremented_prod_index<Seq,range_to_decrement>::decremented_Seq decremented_Seq;
    typedef typename boost::mpl::insert_range<erased_Seq,typename boost::mpl::end<erased_Seq>::type,decremented_Seq>::type newSeq;


};

template<class Seq>
struct decrement_back_indices<Seq,0>
{
    typedef Seq newSeq;
};

}

template<class l_Seq,class r_Seq>
struct MIAPullIndicesUtil
{
    typedef typename internal::pull_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value>::outer_product_indices l_indices;
    typedef typename internal::pull_product_indices<r_Seq,l_Seq,boost::mpl::empty<r_Seq>::value>::outer_product_indices r_indices;
    typedef typename internal::pull_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value>::inter_product_indices inter_product_indices;
    typedef typename boost::mpl::insert_range<l_indices,typename boost::mpl::end<l_indices>::type,r_indices>::type concat;

    typedef typename boost::mpl::insert_range<concat,typename boost::mpl::end<concat>::type,inter_product_indices>::type final_sequence;
    static constexpr size_t MIA_return_order= boost::mpl::size<final_sequence>::value;
    static constexpr size_t inter_product_number=boost::mpl::size<inter_product_indices>::value;
};

template<class l_Seq,class r_Seq>
struct isPureOuterInterProduct
{
    typedef typename internal::pull_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value>::outer_product_indices outer_l_indices;
    typedef typename internal::pull_product_indices<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value>::inter_product_indices inter_l_indices;
    typedef typename boost::mpl::equal_to<
                typename boost::mpl::plus<
                    typename boost::mpl::size<outer_l_indices>::type,
                    typename boost::mpl::size<inter_l_indices>::type
                >::type,
                typename boost::mpl::size<l_Seq>::type
            >::type type;

};

template<class L_MIA, class R_MIA,class l_Seq,class r_Seq >
struct MIAProductUtil: MIAPullIndicesUtil<l_Seq,r_Seq>
{



    typedef typename boost::mpl::if_<
                    typename isPureOuterInterProduct<l_Seq,r_Seq>::type,
                    typename MIANoLatticeProductReturnType<L_MIA,R_MIA,MIAPullIndicesUtil<l_Seq,r_Seq>::MIA_return_order>::type,
                    typename MIAProductReturnType<L_MIA,R_MIA,MIAPullIndicesUtil<l_Seq,r_Seq>::MIA_return_order>::type
                 >::type MIA_return_type;


};

template<class L_MIA, class R_MIA,class l_Seq,class r_Seq >
struct MIASolveUtil: MIAPullIndicesUtil<l_Seq,r_Seq>
{




    typedef typename MIASolveReturnType<L_MIA,R_MIA,MIAPullIndicesUtil<l_Seq,r_Seq>::MIA_return_order>::type MIA_return_type;

};




/*! @} */
/*! @} */

}


#endif // EXPRUTIL_H
