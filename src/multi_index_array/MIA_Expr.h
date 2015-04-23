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
#include <boost/mpl/equal.hpp>
#include <boost/preprocessor/repetition/enum.hpp>

#include "kennytm/vtmp.hpp"
#include "LibMIAUtil.h"
#include "ExprUtil.h"

namespace LibMIA
{







namespace internal
{

/* beginning of taken from KennyTM's answer to this question:
http://stackoverflow.com/questions/10830406/copy-an-mplvector-c-to-a-static-array-at-compile-time*/
template <typename MPLVectorType>
class to_std_array
{
    
	typedef boost::mpl::vector_c<int> numbers;
	//convert to a mpl::vector_c, as the intel compiler seemed to pass mpl::vector instead to this struct, so just make sure we are actually working with an mpl vector_c
	typedef boost::mpl::copy<
		MPLVectorType
		, boost::mpl::back_inserter< numbers >
	>::type result;
	
	static constexpr size_t length = boost::mpl::size<result>::value;
    typedef std::array<int, length> array_type;

    template <size_t... indices>
    static constexpr array_type
            make(const utils::vtmp::integers<indices...>&) noexcept
    {
        return array_type{{
				boost::mpl::at_c<result, indices>::type::value...
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



//template <typename MPLVectorType>
//class to_std_array
//{
//	
//	static constexpr size_t length = boost::mpl::size<MPLVectorType>::value;
//	typedef std::array<int, length> array_type;
//
//	
//
//public:
//	static constexpr array_type make() noexcept
//	{
//		#define LIBMIA_ARRAY_CONVERSION_MACRO(z, i, data) boost::mpl::at_c<data,i>::value
//		static const array_type _array = { BOOST_PP_ENUM(N, LIBMIA_ARRAY_CONVERSION_MACRO, MPLVectorType) };
//		return _array;
//	}
//};


}



//helper class to perform cartesian check on two mpl::vector sequences of indices. An optional Pred template parameter allows the specification of different
//equality prediacates (ie, just check that the integer ID in the indices match, and ignore the ElemWise counter. You must also specify how many matches are allowed.
//For instance, if an MIA product is performed, only one match is allowed between two indices specifying an inner or element-wise product.
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

//helper class to perform an auto check on a mpl::vector sequence of indices. An optional Pred template parameter allows the specification of different
//equality prediacates (ie, just check that the integer ID in the indices match, and ignore the ElemWise counter)
template<class Seq, int rule,typename Pred=boost::mpl::quote2<internal::same_product_index_id> >
struct perform_auto_check
{
    typedef internal::auto_cartesian_product_indices<Seq,boost::mpl::empty<Seq>::value,rule,Pred> sequence_check;
    static void run()
    {
        static_assert(sequence_check::value,"Repeated index in operand.");
    }

};


//helper class to perform merge operations (addition or subtraction). If the indices between the two MIAs match completely, the class will delegate to a
//simpler merge operation that doesn't perform any index shuffle calculations
template<class l_Seq,class r_Seq,class Enable = void>
struct perform_merge
{
};
template<class l_Seq,class r_Seq>
struct perform_merge<l_Seq,r_Seq,typename boost::enable_if< boost::mpl::equal<l_Seq,r_Seq,boost::mpl::quote2<internal::same_product_index_id>> >::type>
{


    //if the two sequences match
    template<class lMIAType,class rMIAType,class Op>
    static auto run( lMIAType & lMIA,  rMIAType & rMIA, const Op & op)->
            typename MIAMergeReturnType<lMIAType,rMIAType>::type *
    {
        typedef typename MIAMergeReturnType<lMIAType,rMIAType>::type cType;

        cType* cMIA(new cType(lMIA.outside_merge(rMIA,op)));


        return cMIA;

    }

};
template<class l_Seq,class r_Seq>
struct perform_merge<l_Seq,r_Seq,typename boost::disable_if< boost::mpl::equal<l_Seq,r_Seq,boost::mpl::quote2<internal::same_product_index_id>>>::type>
{
    //if the two sequences don't match
    template<class lMIAType,class rMIAType,class Op>
    static auto run( lMIAType & lMIA,  rMIAType & rMIA,const Op & op)->
            typename MIAMergeReturnType<lMIAType,rMIAType>::type *

    {


        static_assert(internal::order<typename std::remove_const<lMIAType>::type>::value==internal::order<typename std::remove_const<rMIAType>::type >::value,"Orders of two MIAs must be the same to perform addition.");

        //statically verify MIA indices match correctly
        typedef perform_cartesian_check<l_Seq,r_Seq,internal::merge_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated *EDIT - when unary operations are implemented, we can operate under the assumption that any repeated indices are dealt with
        //perform_auto_check<l_Seq,internal::binary_rule>::run();

        typedef typename MIAMergeReturnType<lMIAType,rMIAType>::type cType;
        //from the mpl::vector of indices, create an mpl::vector of integers specifying how B's indices are shuffled based on A's
        typedef internal::pull_match_order<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;



        cType* cMIA(new cType(lMIA.outside_merge(rMIA,op,internal::to_std_array<typename pulling_index_order::match_order>::make())));
        return cMIA;


    }


};


//same as the above structure, except specialized for destructive mergers (e.g., +=)
template<class l_Seq,class r_Seq,class Enable = void>
struct perform_destructive_merge
{
};
template<class l_Seq,class r_Seq>
struct perform_destructive_merge<l_Seq,r_Seq,typename boost::enable_if< boost::mpl::equal<l_Seq,r_Seq,boost::mpl::quote2<internal::same_product_index_id>>>::type>
{

    //if the two sequences match
    template<class lMIAType,class rMIAType,class Op>
    static void run(lMIAType & lMIA,  rMIAType & rMIA,const Op& op)

    {



        lMIA.merge(rMIA,op);



    }

};
template<class l_Seq,class r_Seq>
struct perform_destructive_merge<l_Seq,r_Seq,typename boost::disable_if< boost::mpl::equal<l_Seq,r_Seq,boost::mpl::quote2<internal::same_product_index_id>> >::type>
{
    //if the two sequences don't match
    template<class lMIAType,class rMIAType,class Op>
    static void run(lMIAType & lMIA,  rMIAType & rMIA,const Op& op)
    {


        static_assert(internal::order<typename std::remove_const<lMIAType>::type>::value==internal::order<typename std::remove_const<rMIAType>::type >::value,"Orders of two MIAs must be the same to perform addition.");

        typedef perform_cartesian_check<l_Seq,r_Seq,internal::merge_rule,boost::mpl::quote2<internal::same_product_index_id> > cartesian_check;
        cartesian_check::run();

        //check to makes sure no left-hand indice is repeated
        perform_auto_check<l_Seq,internal::binary_rule>::run();

        typedef internal::pull_match_order<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value,boost::mpl::quote2<internal::same_product_index_id>> pulling_index_order;



        lMIA.merge(rMIA,op,internal::to_std_array<typename pulling_index_order::match_order>::make());


    }


};

//helper structure to perform index decrements. ie, turn an index sequence (i,!j,!!k) to (i,j,!k). In LibMIA, this operation is performed using ~a(i,!j,!!k).
//the inter_product_number parameter, specifies where the indices that need to be decermented start (from the back). This is valid, because the ~ operation should
//only be performed after an MIA multiplication or solve, which places any interproduct indices at the back after the operation. If inter_product_number is zero,
//this structure delegates to a method that does nothing
struct index_decrementer
{

    //decrement indices method
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

    //do nothing method
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


//delegates to the MIA's toLatticeDiscard or toLatticeExpression method. They are used when the MIA is a temporary or permanent value respectively,
//the former occurs if the MIA is a temporary created as part of an expression. In this case the MIA_Atom has ownership of the MIA. For DenseMIAs,
//the lattice creation method can use an in-place permute instead of a copy
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


//main job of this class is to either call lattice making functions that require a permute or if the MIA data is already arrayed correctly for the lattice,
//then just lattice making method that wraps MIA data directly (choice performed at compile time, hence its ugliness)
struct lattice_permutation_delegator
{

    //if the left MIA is mapped to a lattice using a non-consecutive set of indices, we need to permute
    template<typename _MIA, typename helper,bool hasOwnership,
        typename boost::disable_if<typename helper::is_left_consecutive,int>::type=0
    >
    static auto left_lattice_apply(_MIA & _mia)
    ->decltype(lattice_maker<_MIA,decltype(helper::left_row_order()),decltype(helper::left_column_order()),decltype(helper::left_tab_order()),hasOwnership>
            ::apply(_mia,helper::left_row_order(),helper::left_column_order(),helper::left_tab_order()))
    {


        return lattice_maker<_MIA,decltype(helper::left_row_order()),decltype(helper::left_column_order()),decltype(helper::left_tab_order()),hasOwnership>
            ::apply(_mia,helper::left_row_order(),helper::left_column_order(),helper::left_tab_order());
    }

    //if the left MIA is mapped to a lattice using a consecutive set of indices, no permute
    template<typename _MIA, typename helper,bool hasOwnership,
        typename boost::enable_if<typename helper::is_left_consecutive,int>::type=0
    >
    static auto left_lattice_apply(_MIA & _mia)
    ->decltype(_mia.toStraightLattice(  boost::mpl::size<
                                            typename helper::left_row_inds
                                        >::value,
                                        boost::mpl::size<
                                            typename helper::left_column_inds
                                        >::value
                                    )
               )
    {



        return _mia.toStraightLattice(  boost::mpl::size<
                                            typename helper::left_row_inds
                                        >::value,
                                        boost::mpl::size<
                                            typename helper::left_column_inds
                                        >::value
                                        );
    }
    //right lattice

    //if the right MIA is mapped to a lattice using a non-consecutive set of indices, we need to permute
    template<typename _MIA, typename helper,bool hasOwnership,
        typename boost::disable_if<typename helper::is_right_consecutive,int>::type=0
    >
    static auto right_lattice_apply(_MIA & _mia)
    ->decltype(lattice_maker<_MIA,decltype(helper::right_row_order()),decltype(helper::right_column_order()),decltype(helper::right_tab_order()),hasOwnership>
            ::apply(_mia,helper::right_row_order(),helper::right_column_order(),helper::right_tab_order()))
    {


        return lattice_maker<_MIA,decltype(helper::right_row_order()),decltype(helper::right_column_order()),decltype(helper::right_tab_order()),hasOwnership>
            ::apply(_mia,helper::right_row_order(),helper::right_column_order(),helper::right_tab_order());
    }

    //if the right MIA is mapped to a lattice using a consecutive set of indices, no permute
    template<typename _MIA, typename helper,bool hasOwnership,
        typename boost::enable_if<typename helper::is_right_consecutive,int>::type=0
    >
    static auto right_lattice_apply(_MIA & _mia)
    ->decltype(_mia.toStraightLattice(  boost::mpl::size<
                                            typename helper::right_row_inds
                                        >::value,
                                        boost::mpl::size<
                                            typename helper::right_column_inds
                                        >::value
                                    )
               )
    {
        return _mia.toStraightLattice(  boost::mpl::size<
                                            typename helper::right_row_inds
                                        >::value,
                                        boost::mpl::size<
                                            typename helper::right_column_inds
                                        >::value
                                        );
    }
};







//given three boost::mpl::vectors of integer sequences, checks whether they express a consecutive sequences of integers starting from 0.
template<class first,class second,class third>
struct consecutive_sequence_checker{



    //the total consequetive order of outer product, inner product, and inter product indices
    typedef typename boost::mpl::insert_range<first,typename boost::mpl::end<first>::type, second>::type first_concat;
    typedef typename boost::mpl::insert_range<first_concat,typename boost::mpl::end<first_concat>::type,third>::type final_sequence;
    //is the order of the lattice indices consequetive? If so, we can skip permutation of data elements
    typedef typename boost::mpl::equal<
        boost::mpl::range_c<
            int,0,boost::mpl::size<final_sequence>::value
        >,
        final_sequence,
		boost::mpl::equal_to<boost::mpl::placeholders::_1, boost::mpl::placeholders::_2>
    >::type is_consecutive;

};


//give two boost::mpl::vectors of ProdInds, pulls the integer order of how each sequence matches with the other
template<class l_Seq,class r_Seq>
struct solve_product_expr_helper
{


    typedef typename internal::pull_left_operand_index_sequence<l_Seq,r_Seq,boost::mpl::empty<l_Seq>::value,0> left_index_order;
    typedef typename internal::pull_left_operand_index_sequence<r_Seq,l_Seq,boost::mpl::empty<r_Seq>::value,0> right_outer_index_order;
    typedef typename internal::pull_right_index_order<l_Seq,r_Seq,boost::mpl::empty<r_Seq>::value> pulling_index_order;

    //Assume the following two indices are given l_Seq=(!i,j,k,!n,l) r_Seq=(l,!n,j,m,!i)

    typedef typename left_index_order::no_match_order_sequence left_outer_product_idx; //consecutive list of outer product indices of the left MIA {2}
    typedef typename left_index_order::match_order_sequence left_inner_product_idx; //consecutive list of inner product indices of the left MIA {1,4}
    typedef typename left_index_order::inter_match_order_sequence left_inter_product_idx; //consecutive list of inter product indices of the left MIA {0,3}

    typedef typename right_outer_index_order::no_match_order_sequence right_outer_product_idx; //consecutive list of outer product indices of the right MIA {3}
    typedef typename pulling_index_order::match_order right_inner_product_idx; //the order of right's inner product indices when matched with left's consecutive list of inner product indices {2,0}
    typedef typename pulling_index_order::inter_match_order right_inter_product_idx; //the order of right's inter product indices when matched with left's consecutive list of inter product indices {4,1}

    //return the runtime dimensions of the new MIA
    template<class l_MIA_type, class r_MIA_type>
    static std::array<typename l_MIA_type::index_type,
                boost::mpl::size<left_outer_product_idx>::value +
                boost::mpl::size<right_outer_product_idx>::value+
                boost::mpl::size<left_inter_product_idx>::value>
        run(const l_MIA_type& l_MIA, const r_MIA_type& r_MIA)
    {
        //check to makes sure no left-hand indice is repeated
        perform_auto_check<l_Seq,internal::binary_rule>::run();
        //check to makes sure no right-hand indice is repeated
        perform_auto_check<r_Seq,internal::binary_rule>::run();
        perform_cartesian_check<l_Seq,r_Seq,internal::product_rule>::run();
        //make sure indices that share the same id don't have different elemwise counters
        perform_cartesian_check<l_Seq,r_Seq,internal::product_rule_different_elem_wise,boost::mpl::quote2<internal::same_product_index_id_different_elem_counters>>::run();



        //        print_array(left_outer_product_order, "left_outer");
        //        print_array(left_inner_product_order,"left_inner");
        //        print_array(left_inter_product_order,"left_inter");
        //
        //        print_array(right_outer_product_order, "right_outer");
        //        print_array(right_inner_product_order,"right_inner");
        //        print_array(right_inter_product_order,"right_inter");
        std::array<typename l_MIA_type::index_type,
               boost::mpl::size<left_outer_product_idx>::value +
                boost::mpl::size<right_outer_product_idx>::value+
                boost::mpl::size<left_inter_product_idx>::value> cMIA_dims;
        size_t curIdx=0;
        internal::reorder_from(l_MIA.dims(),left_outer_product_order(),cMIA_dims,curIdx);
        internal::reorder_from(r_MIA.dims(),right_outer_product_order(),cMIA_dims,curIdx);
        internal::reorder_from(l_MIA.dims(),left_inter_product_order(),cMIA_dims,curIdx);
        return cMIA_dims;

    }
    //following methods create a runtime std::array object from the static boost::mpl::vectors of index matching orders
    static constexpr auto left_outer_product_order()->decltype(internal::to_std_array<left_outer_product_idx>::make())
    {
        return internal::to_std_array<left_outer_product_idx>::make();
    }
    static constexpr auto right_outer_product_order()->decltype(internal::to_std_array<right_outer_product_idx>::make())
    {
        return internal::to_std_array<right_outer_product_idx>::make();
    }
    static constexpr auto left_inner_product_order()->decltype(internal::to_std_array<left_inner_product_idx>::make())
    {
        return internal::to_std_array<left_inner_product_idx>::make();
    }
    static constexpr auto right_inner_product_order()->decltype(internal::to_std_array<right_inner_product_idx>::make())
    {
        return internal::to_std_array<right_inner_product_idx>::make();
    }
    static constexpr auto left_inter_product_order()->decltype(internal::to_std_array<left_inter_product_idx>::make())
    {
        return internal::to_std_array<left_inter_product_idx>::make();
    }
    static constexpr auto right_inter_product_order()->decltype(internal::to_std_array<right_inter_product_idx>::make())
    {
        return internal::to_std_array<right_inter_product_idx>::make();
    }
};

//chooses which indices map to rows/columns/tabs (for solving left outerproduct are mapped to columns and inner product to rows)
template<class mia_expr_helper>
struct solve_lattice_expr_helper
{
    typedef typename mia_expr_helper::left_inner_product_idx left_row_inds;
    typedef typename mia_expr_helper::left_outer_product_idx left_column_inds;
    typedef typename mia_expr_helper::left_inter_product_idx left_tab_inds;
    typedef typename consecutive_sequence_checker<left_row_inds,left_column_inds,left_tab_inds>::is_consecutive is_left_consecutive;

    typedef typename mia_expr_helper::right_inner_product_idx right_row_inds;
    typedef typename mia_expr_helper::right_outer_product_idx right_column_inds;
    typedef typename mia_expr_helper::right_inter_product_idx right_tab_inds;
    typedef typename consecutive_sequence_checker<right_row_inds,right_column_inds,right_tab_inds>::is_consecutive is_right_consecutive;

};

//chooses which indices map to rows/columns/tabs (for solving left outerproduct are mapped to rows and inner product to columns)
template<class mia_expr_helper>
struct product_lattice_expr_helper
{
    typedef typename mia_expr_helper::left_outer_product_idx left_row_inds;
    typedef typename mia_expr_helper::left_inner_product_idx left_column_inds;
    typedef typename mia_expr_helper::left_inter_product_idx left_tab_inds;
    typedef typename consecutive_sequence_checker<left_row_inds,left_column_inds,left_tab_inds>::is_consecutive is_left_consecutive;

    typedef typename mia_expr_helper::right_inner_product_idx right_row_inds;
    typedef typename mia_expr_helper::right_outer_product_idx right_column_inds;
    typedef typename mia_expr_helper::right_inter_product_idx right_tab_inds;
    typedef typename consecutive_sequence_checker<right_row_inds,right_column_inds,right_tab_inds>::is_consecutive is_right_consecutive;

};



//helper class to return the row, column, and tab index orders (similar to what solve_product_expr_helper does for inner, outer, and inter products).
//However, the class is paramterized by super_lattice_expr_helper, which delegates indices to rows or columns differently depending on whether a solve
//or product is being performed.
//superclass should be product_lattice_expr_helper or solve_lattice_expr_helper
template<class super_lattice_expr_helper>
struct lattice_expr_helper : public super_lattice_expr_helper
{

    static constexpr auto left_row_order()->decltype(internal::to_std_array<typename super_lattice_expr_helper::left_row_inds>::make())
    {
        return internal::to_std_array<typename super_lattice_expr_helper::left_row_inds>::make();
    }
    static constexpr auto right_row_order()->decltype(internal::to_std_array<typename super_lattice_expr_helper::right_row_inds>::make())
    {
        return internal::to_std_array<typename super_lattice_expr_helper::right_row_inds>::make();
    }
    static constexpr auto left_column_order()->decltype(internal::to_std_array<typename super_lattice_expr_helper::left_column_inds>::make())
    {
        return internal::to_std_array<typename super_lattice_expr_helper::left_column_inds>::make();
    }
    static constexpr auto right_column_order()->decltype(internal::to_std_array<typename super_lattice_expr_helper::right_column_inds>::make())
    {
        return internal::to_std_array<typename super_lattice_expr_helper::right_column_inds>::make();
    }
    static constexpr auto left_tab_order()->decltype(internal::to_std_array<typename super_lattice_expr_helper::left_tab_inds>::make())
    {
        return internal::to_std_array<typename super_lattice_expr_helper::left_tab_inds>::make();
    }
    static constexpr auto right_tab_order()->decltype(internal::to_std_array<typename super_lattice_expr_helper::right_tab_inds>::make())
    {
        return internal::to_std_array<typename super_lattice_expr_helper::right_tab_inds>::make();
    }


};


template<class Seq>
struct contraction_helper
{
    typedef internal::pull_contract_indices<Seq,boost::mpl::empty<Seq>::value> contraction_index_helper;
    typedef internal::pull_sequence_from_list<Seq,typename contraction_index_helper::contract_indices,boost::mpl::empty<typename contraction_index_helper::contract_indices>::value> index_sequences;
    typedef typename contraction_index_helper::remaining_indices remaining_indices;
    typedef typename boost::mpl::size<remaining_indices>::type remaining_indices_size;
    typedef typename contraction_index_helper::contract_indices contract_indices;
    typedef typename boost::mpl::size<contract_indices>::type contract_indices_size;
    typedef typename index_sequences::sequence contraction_index_sequence;
    typedef typename boost::mpl::size<contraction_index_sequence>::type contraction_index_sequence_size;
    typedef typename index_sequences::partition_sequence contraction_partition_sequence;

};

template<class Seq>
struct attraction_helper
{
    typedef internal::pull_attract_indices<Seq,boost::mpl::empty<Seq>::value> attraction_index_helper;
    typedef internal::pull_sequence_from_list<Seq,typename attraction_index_helper::attract_indices,boost::mpl::empty<typename attraction_index_helper::attract_indices>::value> index_sequences;

    typedef typename attraction_index_helper::attract_indices attract_indices;
    typedef typename boost::mpl::size<attract_indices>::type attract_indices_size;
    typedef typename index_sequences::sequence attraction_index_sequence;
    typedef typename boost::mpl::size<attraction_index_sequence>::type attraction_index_sequence_size;
    typedef typename index_sequences::partition_sequence attraction_partition_sequence;

};

template<class Seq>
struct unary_helper
{

    typedef internal::pull_contract_indices<Seq,boost::mpl::empty<Seq>::value> contraction_index_helper;
    typedef typename contraction_index_helper::remaining_indices remaining_indices_after_contraction; //mpl::vector of MIA indices not undergoing contraction
    //get mpl::vector of indices undergoing contraction
    typedef internal::pull_sequence_from_list<Seq,typename contraction_index_helper::contract_indices,boost::mpl::empty<typename contraction_index_helper::contract_indices>::value> contract_index_puller;
    typedef typename contract_index_puller::sequence contraction_index_sequence;
    typedef typename boost::mpl::size<contraction_index_sequence>::type contraction_index_sequence_size;
    //get mpl::vector of partitions of the possibly more than one contraction set
    typedef typename contract_index_puller::partition_sequence contraction_partition_sequence;
    //pull the attraction indices from the indices remaining after contraction
    typedef internal::pull_attract_indices<remaining_indices_after_contraction,boost::mpl::empty<remaining_indices_after_contraction>::value> attraction_index_helper;
    typedef typename attraction_index_helper::remaining_indices remaining_indices_after_attraction_and_contraction; //mpl::vector of MIA indices not undergoing contraction or attraction

    //get the location of these attraction indices and their partitions
    typedef internal::pull_sequence_from_list<Seq,typename attraction_index_helper::attract_indices,boost::mpl::empty<typename attraction_index_helper::attract_indices>::value> attract_index_puller;
    typedef typename attract_index_puller::sequence attraction_index_sequence;
    typedef typename boost::mpl::size<attraction_index_sequence>::type attraction_index_sequence_size;
    typedef typename attract_index_puller::partition_sequence attraction_partition_sequence;

    //add the attraction indices to the end of the list of indices remaining after contraction to get a combined list of remaining indices
    typedef typename boost::mpl::insert_range<remaining_indices_after_attraction_and_contraction,
                                                    typename boost::mpl::end<remaining_indices_after_attraction_and_contraction>::type,
                                                    typename attraction_index_helper::attract_indices
                                                >::type remaining_indices;
    typedef typename boost::mpl::size<remaining_indices>::type remaining_indices_size;

};


//!Only enabled when no contraction or attraction is taking place (known by examining Seq template parameter). Will just return the given MIA with the same Seq, doing nothing else.
template<class _MIA,class Seq,typename boost::disable_if<
                                boost::mpl::less<
                                    typename unary_helper<Seq>::remaining_indices_size,
                                    typename boost::mpl::size<Seq>::type
                                >,
                                int
                            >::type=0
>
MIA_Atom<
    _MIA,
    Seq,
    false
>
perform_unary(_MIA & mia)
{



    return MIA_Atom<
                _MIA,
                Seq,
                false
            >(&mia);


}



//!Only enabled when a contraction or attraction is actually taking place (known by examining Seq template parameter). Will return a new MIA (of smaller order) with a new index sequence.
template<class _MIA,class Seq,typename boost::enable_if<
                                boost::mpl::less<
                                    typename unary_helper<Seq>::remaining_indices_size,
                                    typename boost::mpl::size<Seq>::type
                                >,
                                int
                            >::type=0
>
MIA_Atom<
    typename MIAUnaryType<_MIA,unary_helper<Seq>::remaining_indices_size::value>::type,
    typename unary_helper<Seq>::remaining_indices,
    true
>
perform_unary(_MIA & mia)
{


//    print_array(internal::to_std_array<typename unary_helper<Seq>::contraction_index_sequence>::make(),"contraction_index_sequence");
//    print_array(internal::to_std_array<typename unary_helper<Seq>::contraction_partition_sequence>::make(),"contraction_partition_sequence");
//    print_array(internal::to_std_array<typename unary_helper<Seq>::attraction_index_sequence>::make(),"attraction_index_sequence");
//    print_array(internal::to_std_array<typename unary_helper<Seq>::attraction_partition_sequence>::make(),"attraction_partition_sequence");


    //print_array(internal::to_std_array<typename contraction_index_helper::contract_sequence>::make(),"contract indices");

    typedef typename MIAUnaryType<_MIA,unary_helper<Seq>::remaining_indices_size::value>::type retType;
    retType * ret(new retType());


    *ret=mia.contract_attract(internal::to_std_array<typename contraction_helper<Seq>::contraction_index_sequence>::make(),
                              internal::to_std_array<typename contraction_helper<Seq>::contraction_partition_sequence>::make(),
                              internal::to_std_array<typename attraction_helper<Seq>::attraction_index_sequence>::make(),
                              internal::to_std_array<typename attraction_helper<Seq>::attraction_partition_sequence>::make());
    return MIA_Atom<
                retType,
                typename unary_helper<Seq>::remaining_indices,
                true
            >(ret);


}


//!Expression class created whenever an MIA is index by ProdInds. Reponsible for all compile-time delegation of operations
template<class _MIA,class m_Seq,bool ownership,size_t inter_product_number>
class MIA_Atom
{
public:
    _MIA* m_mia;
    static constexpr bool mHasOwnership=ownership;
    MIA_Atom(_MIA* mia): m_mia(mia)
    {
        static_assert(internal::is_MIA<typename std::remove_const<_MIA>::type>::value,"Somehow expression was instantiated with a non-MIA class.");


    }
    ~MIA_Atom(){
        if(m_mia && mHasOwnership)
            delete m_mia;

    }


    operator _MIA &&(){
        return std::move(*(this->m_mia));
    }

    //! When certain operations (such as multiplication) are performed using two of the same MIAs, a copy needs to be done
    /*! Only need to even consider making a copy if the MIA_Atom doesn't own the MIA (otherwise it's a temporary) and if
        it's a sparse MIA, because SparseMIAs are the only kind that will alter their data (through a sort) even when the
        MIA_Atom doesn't own the MIA. An MIA Atom should not be able to accept a const SparseMIA, so the cast to void *
        (instead of const void *) should be fine.
    */
    template<bool fromOwnership,class fromMIAType, class checkMIAType,
        typename boost::enable_if<
            boost::mpl::and_<
                boost::mpl::bool_<!fromOwnership>,
                boost::mpl::and_<
                    internal::is_SparseMIA<typename std::remove_const<fromMIAType>::type>,
                    internal::is_SparseMIA<typename std::remove_const<checkMIAType>::type>
                >
            >,
            int
        >::type=0
    >
    bool check_if_copy_needed(fromMIAType * _from, checkMIAType * _check){
        if(static_cast<void*>(_from)==static_cast<void*>(_check))
            return true;
        return false;

    }

    template<bool fromOwnership,class fromMIAType, class checkMIAType,
        typename boost::disable_if<
            boost::mpl::and_<
                boost::mpl::bool_<!fromOwnership>,
                boost::mpl::and_<
                    internal::is_SparseMIA<typename std::remove_const<fromMIAType>::type>,
                    internal::is_SparseMIA<typename std::remove_const<checkMIAType>::type>
                >
            >,
            int
        >::type=0
    >
    bool check_if_copy_needed(fromMIAType * _from, checkMIAType * _check){

        return false;

    }

    //should only be enabled when we've reduce two MIAs to a single number (ie pure inner product), in this case it just returns a POD number
    template<class otherMIA,class r_Seq,bool otherOwnership,size_t other_inter_number,
                typename boost::enable_if_c<
                    boost::mpl::size<
                        typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence
                    >::value==0,
                    int
                >::type=0
            >
    auto operator*(const MIA_Atom<otherMIA,r_Seq, otherOwnership,other_inter_number> & Rhs)->
        typename internal::data_type<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type
        >::type
    {


        typedef solve_product_expr_helper<m_Seq,r_Seq> mia_expr_helper;
#ifdef LIBMIA_CHECK_DIMS

        check_dims_match(m_mia->dims(),mia_expr_helper::left_inner_product_order(),Rhs.m_mia->dims(),mia_expr_helper::right_inner_product_order(),"An inner-product index is matched with a different sized index.");
        check_dims_match(m_mia->dims(),mia_expr_helper::left_inter_product_order(),Rhs.m_mia->dims(),mia_expr_helper::right_inter_product_order(),"An inter-product index is matched with a different sized index.");

#endif

        typedef lattice_expr_helper<product_lattice_expr_helper<mia_expr_helper>> m_lattice_expr_helper;

        auto aLat=lattice_permutation_delegator::left_lattice_apply<_MIA,m_lattice_expr_helper,mHasOwnership>(*m_mia);
        //std::cout << "ALat done" << std::endl;
        //aLat.print();

        //std::cout << "BLat done" << std::endl;
        //bLat.print();



        if(check_if_copy_needed<otherOwnership>(m_mia,Rhs.m_mia)){ //is the mia type is sparse and they refer to the same object, then we must make a copy, because the lattice making process uses sort rather than a copy
            typedef typename MIANonlinearFuncType<otherMIA>::type bMIACopyType; //if the MIA is a 'mapped' datatype, we need to copy it to a MIA type that owns its data
            bMIACopyType temp(*(Rhs.m_mia));
            auto bLat=lattice_permutation_delegator::right_lattice_apply<bMIACopyType,m_lattice_expr_helper,true>(temp);
            auto cLat=aLat*bLat;

            return *(cLat.data_begin());
        }
        else{
            auto bLat=lattice_permutation_delegator::right_lattice_apply<otherMIA,m_lattice_expr_helper,otherOwnership>(*(Rhs.m_mia));
            auto cLat=aLat*bLat;
            return *(cLat.data_begin());
        }


//        std::cout << "CLat done" << std::endl;
//        cLat.print();




    }



    //when performing a lattice product
    template<class otherMIA,class r_Seq,bool otherOwnership,size_t other_inter_number,
        typename boost::disable_if<
            boost::mpl::or_<
                typename isPureOuterInterProduct<m_Seq,r_Seq>::type,
                boost::mpl::bool_<
                    boost::mpl::size<
                        typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence
                    >::value==0
                >
            >,
            int
        >::type=0
    >
    auto operator*(const MIA_Atom<otherMIA,r_Seq, otherOwnership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            true,
            MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    {

        typedef solve_product_expr_helper<m_Seq,r_Seq> mia_expr_helper;
#ifdef LIBMIA_CHECK_DIMS

        check_dims_match(m_mia->dims(),mia_expr_helper::left_inner_product_order(),Rhs.m_mia->dims(),mia_expr_helper::right_inner_product_order(),"An inner-product index is matched with a different sized index.");
        check_dims_match(m_mia->dims(),mia_expr_helper::left_inter_product_order(),Rhs.m_mia->dims(),mia_expr_helper::right_inter_product_order(),"An inter-product index is matched with a different sized index.");

#endif
        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;

        auto cMIA_dims=mia_expr_helper::run(*m_mia,*(Rhs.m_mia));

        typedef lattice_expr_helper<product_lattice_expr_helper<mia_expr_helper>> m_lattice_expr_helper;


        auto aLat=lattice_permutation_delegator::left_lattice_apply<_MIA,m_lattice_expr_helper,mHasOwnership>(*m_mia);

        constexpr size_t _inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
        MIA_return_type* cMIA;
        if(check_if_copy_needed<otherOwnership>(m_mia,Rhs.m_mia)){ //is the mia type is sparse and they refer to the same object, then we must make a copy, because the lattice making process uses sort rather than a copy
            typedef typename MIANonlinearFuncType<otherMIA>::type bMIACopyType; //if the MIA is a 'mapped' datatype, we need to copy it to a MIA type that owns its data
            bMIACopyType temp(*(Rhs.m_mia));
            auto bLat=lattice_permutation_delegator::right_lattice_apply<bMIACopyType,m_lattice_expr_helper,true>(temp);
            auto cLat=aLat*bLat;

            cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
        }
        else{
            auto bLat=lattice_permutation_delegator::right_lattice_apply<otherMIA,m_lattice_expr_helper,otherOwnership>(*(Rhs.m_mia));
            auto cLat=aLat*bLat;
            cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
        }


        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,_inter_product_number>(cMIA);




    }


//    //!Used to actually perform lattice products
//    /*!Normally we can just perform a normal lattice product. But if we have two sparse MIAs, we can delegate to a poly-algorithm that
//    performs different sparse matrix computations depending on the hyper-sparsity of the operands. This function is called when we do
//    NOT have two sparse MIAs
//
//    */
//    template<class otherMIA,class r_Seq,bool otherOwnership,size_t other_inter_number,
//        typename boost::disable_if<
//            boost::mpl::and_<
//                internal::is_SparseMIA<typename std::remove_const<_MIA>::type>,
//                internal::is_SparseMIA<typename std::remove_const<otherMIA>::type>
//            >,
//            int
//        >::type=0
//    >
//    auto doProduct(otherMIA & r_mia, MIA_return_type* cMIA )->
//        typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type *
//    {
//        typedef solve_product_expr_helper<m_Seq,r_Seq> mia_expr_helper;
//        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
//
//        auto cMIA_dims=mia_expr_helper::run(*m_mia,*(Rhs.m_mia));
//
//        typedef lattice_expr_helper<product_lattice_expr_helper<mia_expr_helper>> m_lattice_expr_helper;
//
//
//        auto aLat=lattice_permutation_delegator::left_lattice_apply<_MIA,m_lattice_expr_helper,mHasOwnership>(*m_mia);
//
//        constexpr size_t _inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
//        MIA_return_type* cMIA;
//
//        auto bLat=lattice_permutation_delegator::right_lattice_apply<otherMIA,m_lattice_expr_helper,otherOwnership>(*(Rhs.m_mia));
//        auto cLat=aLat*bLat;
//        cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
//
//        return cMIA;
//    }
//
//    //!Used to actually perform lattice products
//    /*!Normally we can just perform a normal lattice product. But if we have two sparse MIAs, we can delegate to a poly-algorithm that
//    performs different sparse matrix computations depending on the hyper-sparsity of the operands. This function is called when we have two
//    sparse MIAs
//
//    */
//    template<class otherMIA,class r_Seq,bool otherOwnership,size_t other_inter_number,
//        typename boost::enable_if<
//            boost::mpl::and_<
//                internal::is_SparseMIA<typename std::remove_const<_MIA>::type>,
//                internal::is_SparseMIA<typename std::remove_const<otherMIA>::type>
//            >,
//            int
//        >::type=0
//    >
//    auto doProduct(otherMIA & r_mia, MIA_return_type* cMIA )->
//        typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type *
//    {
//        typedef solve_product_expr_helper<m_Seq,r_Seq> mia_expr_helper;
//        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
//
//        auto cMIA_dims=mia_expr_helper::run(*m_mia,*(Rhs.m_mia));
//
//        typedef lattice_expr_helper<product_lattice_expr_helper<mia_expr_helper>> m_lattice_expr_helper;
//
//        constexpr size_t _inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;
//        MIA_return_type* cMIA;
//        if(check_if_copy_needed<otherOwnership>(m_mia,Rhs.m_mia)){ //is the mia type is sparse and they refer to the same object, then we must make a copy, because the lattice making process uses sort rather than a copy
//            typedef typename MIANonlinearFuncType<otherMIA>::type bMIACopyType; //if the MIA is a 'mapped' datatype, we need to copy it to a MIA type that owns its data
//            bMIACopyType temp(*(Rhs.m_mia));
//            auto cLat=sparseMIAMultPolyAlg(*m_mia,temp,m_lattice_expr_helper::left_row_inds,m_lattice_expr_helper::left_col_inds,m_lattice_expr_helper::left_tab_inds
//                                           m_lattice_expr_helper::right_row_inds,m_lattice_expr_helper::right_col_inds,m_lattice_expr_helper::right_tab_inds);
//
//            cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
//        }
//        else{
//            auto bLat=lattice_permutation_delegator::right_lattice_apply<otherMIA,m_lattice_expr_helper,otherOwnership>(*(Rhs.m_mia));
//            auto cLat=aLat*bLat;
//            cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
//        }
//
//
//        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,_inter_product_number>(cMIA);
//    }


    //when multiplication only has inter or outer products, we don't need to create a lattice, call specialized function in mia
    template<class otherMIA,class r_Seq,bool otherOwnership,size_t other_inter_number,
        typename boost::enable_if<
            boost::mpl::and_<
                typename isPureOuterInterProduct<m_Seq,r_Seq>::type,
                boost::mpl::bool_<
                    boost::mpl::size<
                        typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence
                    >::value!=0
                >
            >,
            int
        >::type=0
    >
    auto operator*(const MIA_Atom<otherMIA,r_Seq, otherOwnership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type,
            typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,
            true,
            MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number
        >
    {


        //std::cout << "Pure outer/inter started " << std::endl;
        typedef typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        typedef solve_product_expr_helper<m_Seq,r_Seq> helper;

#ifdef LIBMIA_CHECK_DIMS

        check_dims_match(m_mia->dims(),helper::left_inner_product_order(),Rhs.m_mia->dims(),helper::right_inner_product_order(),"An inner-product index is matched with a different sized index.");
        check_dims_match(m_mia->dims(),helper::left_inter_product_order(),Rhs.m_mia->dims(),helper::right_inter_product_order(),"An inter-product index is matched with a different sized index.");

#endif

        MIA_return_type* cMIA(new MIA_return_type(m_mia->noLatticeMult(*Rhs.m_mia,helper::left_inter_product_order(),
                                                                       helper::left_outer_product_order(),
                                                                       helper::right_inter_product_order(),
                                                                       helper::right_outer_product_order())));

        constexpr size_t _inter_product_number=MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;

        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIAProductUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,_inter_product_number>(cMIA);

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
        //std::cout << "ALat done" << std::endl;
        //aLat.print();



        typedef solve_product_expr_helper<m_Seq,r_Seq> mia_expr_helper;
        typedef typename MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::MIA_return_type MIA_return_type;
        auto cMIA_dims=mia_expr_helper::run(*m_mia,*(Rhs.m_mia));
        typedef lattice_expr_helper<solve_lattice_expr_helper<mia_expr_helper>> m_lattice_expr_helper;

        auto aLat=lattice_permutation_delegator::left_lattice_apply<_MIA,m_lattice_expr_helper,mHasOwnership>(*m_mia);

        MIA_return_type* cMIA;
        if(check_if_copy_needed<other_ownership>(m_mia,Rhs.m_mia)){ //if the mia type is sparse and they refer to the same object, then we must make a copy, because the lattice making process uses sort rather than a copy
            typedef typename MIANonlinearFuncType<otherMIA>::type bMIACopyType; //if the MIA is a 'mapped' datatype, we need to copy it to a MIA type that owns its data
            bMIACopyType temp(*(Rhs.m_mia));
            auto bLat=lattice_permutation_delegator::right_lattice_apply<bMIACopyType,m_lattice_expr_helper,true>(temp);
            auto cLat=aLat.solve(bLat);

            cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
        }
        else{
            auto bLat=lattice_permutation_delegator::right_lattice_apply<otherMIA,m_lattice_expr_helper,other_ownership>(*(Rhs.m_mia));

            auto cLat=aLat.solve(bLat);
            cMIA=new MIA_return_type(cMIA_dims,std::move(cLat));
        }


        constexpr size_t _inter_product_number=MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::inter_product_number;

        //create an MIA from cLat
        return MIA_Atom<MIA_return_type,typename MIASolveUtil<_MIA,otherMIA,m_Seq,r_Seq>::final_sequence,true,_inter_product_number>(cMIA);




    }




    MIA_Atom& operator=(const MIA_Atom & Rhs)
    {


        //std::cout << "Idential Stuff here " << std::endl;
        if(Rhs.mHasOwnership)
            *m_mia=std::move(*(Rhs.m_mia));
        else
            *m_mia=*(Rhs.m_mia);
        return *this;
    }

    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number,
     typename boost::enable_if< boost::mpl::equal<m_Seq,r_Seq,boost::mpl::quote2<internal::same_product_index_id>>, int >::type = 0
     >
    MIA_Atom& operator=(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)
    {




        static_assert(internal::order<typename std::remove_const<_MIA>::type>::value==internal::order<typename std::remove_const<otherMIA>::type >::value,"Orders of two MIAs must be the same to perform assignment.");

        if(other_ownership){

            *m_mia=std::move(*(Rhs.m_mia));
        }
        else
            *m_mia=*(Rhs.m_mia);
        return *this;



    }

    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number,
     typename boost::disable_if< boost::mpl::equal<m_Seq,r_Seq,boost::mpl::quote2<internal::same_product_index_id>>, int >::type = 0
    >
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



        static_assert(internal::order<typename std::remove_const<_MIA>::type>::value==internal::order<typename std::remove_const<otherMIA>::type>::value,"Orders of two MIAs must be the same to perform addition.");

        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;


//        auto lambda=[](const a_data_type & _a, const  b_data_type & _b){
//            return _a+_b;
//        };

        //perform_destructive_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),lambda);
        perform_destructive_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),std::plus<typename internal::data_type<cType>::type>());





        return *this;


    }


    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number,
        typename boost::disable_if<
            boost::mpl::and_<
                boost::mpl::bool_<mHasOwnership>,
                boost::is_same<
                    _MIA,
                    typename MIAMergeReturnType<_MIA,otherMIA>::type
                >
            >,
            int
        >::type = 0
    >
    auto operator+(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAMergeReturnType<_MIA,otherMIA>::type,
            m_Seq,
            true,
            inter_product_number
        >
    {






        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;

//        auto lambda=[](const a_data_type & _a, const b_data_type & _b){
//            return _a+_b;
//        };

        //cType* cMIA=perform_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),lambda);

        cType* cMIA=perform_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),std::plus<typename internal::data_type<cType>::type>());



        MIA_Atom<cType,m_Seq,true,inter_product_number> C(cMIA);
        return C;

    }

    //if the MIA_Atom owns its MIA and its datatype is the same as the return, then it's temporary, so we can do destructive merge
    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number,
        typename boost::enable_if<
            boost::mpl::and_<
                boost::mpl::bool_<mHasOwnership>,
                boost::is_same<
                    _MIA,
                    typename MIAMergeReturnType<_MIA,otherMIA>::type
                >
            >,
            int
        >::type = 0
    >
    auto operator+(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom
    {



        (*this)+=Rhs;
        MIA_Atom C(this->m_mia);
        this->m_mia=nullptr;
        return C;
    }


    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number>
    MIA_Atom& operator-=(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)
    {



        static_assert(internal::order<typename std::remove_const<_MIA>::type>::value==internal::order<typename std::remove_const<otherMIA>::type>::value,"Orders of two MIAs must be the same to perform subtraction.");

        //check to makes sure no left-hand indice is repeated
        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;

//        auto lambda=[](const a_data_type & _a, const b_data_type & _b){
//            return _a-_b;
//        };

        //perform_destructive_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),lambda);
        perform_destructive_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),std::minus<typename internal::data_type<cType>::type>());






        return *this;


    }


    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number,
        typename boost::disable_if<
            boost::mpl::and_<
                boost::mpl::bool_<mHasOwnership>,
                boost::is_same<
                    _MIA,
                    typename MIAMergeReturnType<_MIA,otherMIA>::type
                >
            >,
            int
        >::type = 0
    >
    auto operator-(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom<
            typename MIAMergeReturnType<_MIA,otherMIA>::type,
            m_Seq,
            true,
            inter_product_number
        >
    {


        typedef typename MIAMergeReturnType<_MIA,otherMIA>::type cType;


//        auto lambda=[](const a_data_type & _a, const b_data_type & _b){
//            return _a-_b;
//        };

        //cType* cMIA=perform_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),lambda);
        cType* cMIA=perform_merge<m_Seq,r_Seq>::template run(*m_mia,*(Rhs.m_mia),std::minus<typename internal::data_type<cType>::type>());




        MIA_Atom<cType,m_Seq,true,inter_product_number> C(cMIA);
        return C;

    }

    //if the MIA_Atom owns its MIA and its datatype is the same as the return, then it's temporary, so we can do destructive merge
    template<class otherMIA,class r_Seq,bool other_ownership,size_t other_inter_number,
        typename boost::enable_if<
            boost::mpl::and_<
                boost::mpl::bool_<mHasOwnership>,
                boost::is_same<
                    _MIA,
                    typename MIAMergeReturnType<_MIA,otherMIA>::type
                >
            >,
            int
        >::type = 0
    >
    auto operator-(const MIA_Atom<otherMIA,r_Seq,other_ownership,other_inter_number> & Rhs)->
        MIA_Atom
    {



        (*this)-=Rhs;
        MIA_Atom C(this->m_mia);
        this->m_mia=nullptr;
        return C;
    }


    auto operator~()->MIA_Atom<_MIA,typename internal::decrement_back_indices<m_Seq,inter_product_number>::newSeq,mHasOwnership>{
        return index_decrementer::template apply<_MIA,m_Seq,mHasOwnership,inter_product_number>(*this);

    }





private:

    //!Helper function to help check whether inner or inter product dimensions match.
    template<size_t size1,size_t size2,size_t shared_size,class index_type1, class index_type2,class index_type3>
    void check_dims_match(const std::array<index_type1,size1> & dims1,const std::array<index_type3,shared_size> & order1,
                     const std::array<index_type2,size2> & dims2,const std::array<index_type3,shared_size> & order2, const std::string & message){

        for(size_t idx=0;idx<shared_size;++idx){
            if (dims1[order1[idx]]!=dims2[order2[idx]])
                throw MIAParameterException(message);
        }

    }


};

//if the MIA_Atom owns the data then its' temporary and we can just modify the MIA data in-place
template<class _MIA,class m_Seq,bool ownership,size_t inter_product_number, class UnaryOperation,typename boost::enable_if_c<ownership,int>::type=0>
MIA_Atom<typename MIANonlinearFuncType<_MIA>::type,m_Seq,true,inter_product_number>
do_func(const MIA_Atom<_MIA,m_Seq,ownership,inter_product_number> & _atom, UnaryOperation op){

    _atom.m_mia->apply_func(op);

    return _atom;
}

//if the MIA_Atom does not the data then we cannot modify data in-place, so just do a copy
template<class _MIA,class m_Seq,bool ownership,size_t inter_product_number, class UnaryOperation,typename boost::disable_if_c<ownership,int>::type=0>
MIA_Atom<typename MIANonlinearFuncType<_MIA>::type,m_Seq,true,inter_product_number>
do_func(const MIA_Atom<_MIA,m_Seq,ownership,inter_product_number> & _atom, UnaryOperation op){
    typedef typename MIANonlinearFuncType<_MIA>::type retType;
    retType * ret=new retType(*(_atom.m_mia));
    ret->apply_func(op);
    return MIA_Atom<retType,m_Seq,true,inter_product_number>(ret);
}



template<class _MIA,class m_Seq,bool ownership,size_t inter_product_number>
MIA_Atom<typename MIANonlinearFuncType<_MIA>::type,m_Seq,true,inter_product_number>
sqrt(const MIA_Atom<_MIA,m_Seq,ownership,inter_product_number> & _atom){
    typedef typename internal::data_type<typename std::remove_const<_MIA>::type>::type data_type;
    data_type (*fpc)(data_type) = &std::sqrt;
    return do_func(_atom,std::ptr_fun(fpc));

}

template<class _MIA,class m_Seq,bool ownership,size_t inter_product_number>
MIA_Atom<typename MIANonlinearFuncType<_MIA>::type,m_Seq,true,inter_product_number>
pow(const MIA_Atom<_MIA,m_Seq,ownership,inter_product_number> & _atom,typename internal::data_type<typename std::remove_const<_MIA>::type>::type  exp){
    typedef typename internal::data_type<typename std::remove_const<_MIA>::type>::type data_type;

    data_type (*fpc)(data_type,data_type) = &std::pow;
    auto f1 = std::bind(fpc, std::placeholders::_1, exp);

    return do_func(_atom,f1);




}

}// namespace LibMIA

#endif // MIA_EXPR_H_INCLUDED
