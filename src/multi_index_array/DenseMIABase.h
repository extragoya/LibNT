// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef DENSEMIABASE_H_INCLUDED
#define DENSEMIABASE_H_INCLUDED

#include <iostream>
#include <array>
#include <algorithm>
#include <numeric>
#include <functional>
#include <boost/dynamic_bitset.hpp>

#include "LibMIAUtil.h"
#include "LibMIARanges.h"
#include "FunctionUtil.h"
#include "MIA.h"
#include "DenseLattice.h"
#include "MappedDenseLattice.h"


namespace LibMIA
{
/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */
namespace internal
{



template<class Derived>
struct data_type<DenseMIABase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct data_type_ref<DenseMIABase<Derived>>: public data_type_ref<Derived> {};

template<class Derived>
struct const_data_type_ref<DenseMIABase<Derived>>: public const_data_type_ref<Derived> {};

template<class Derived>
struct index_type<DenseMIABase<Derived> >: public index_type<Derived> {};

template<class Derived>
struct order<DenseMIABase<Derived> >: public order<Derived> {};

template<class Derived>
struct data_iterator<DenseMIABase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct const_data_iterator<DenseMIABase<Derived> >: public const_data_iterator<Derived> {};



template<class Derived>
struct const_storage_iterator<DenseMIABase<Derived> >: public const_storage_iterator<Derived> {};

template<class Derived>
struct FinalDerived<DenseMIABase<Derived> >:public FinalDerived<Derived>{};

}



//!  Base class for dense multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  dense multi-index array types. Provides operations and functions common to all dense
  multi-index arrays. Due to the existence of ImplictMIA, which is read-only for the raw data, DenseMIABase
  is only meant to define read-only functionality.

  \tparam Derived   should only be a dense multi-index array class.
*/
template <class Derived>
class DenseMIABase: public MIA<DenseMIABase<Derived > >
{
public:


    typedef typename internal::data_type<DenseMIABase>::type data_type;
    typedef typename internal::data_type_ref<DenseMIABase>::type data_type_ref;
    typedef typename internal::const_data_type_ref<DenseMIABase>::type const_data_type_ref;
    typedef typename internal::index_type<DenseMIABase>::type index_type;
    typedef typename internal::data_iterator<DenseMIABase>::type data_iterator;
    typedef typename internal::const_data_iterator<DenseMIABase>::type const_data_iterator;
    typedef typename internal::FinalDerived<DenseMIABase>::type FinalDerived;

    constexpr static size_t mOrder=internal::order<DenseMIABase>::value;
    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }
    /** \returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

    FinalDerived& final_derived()
    {
        return derived().final_derived();
    }
    /** \returns a const reference to the derived object */
    const FinalDerived& final_derived() const
    {
        return derived().final_derived();
    }

    template<typename... Dims>
    DenseMIABase(Dims... dims): MIA<DenseMIABase<Derived > >(dims...) {}


    DenseMIABase(): MIA<DenseMIABase<Derived > >() {}


    DenseMIABase(std::array<index_type,internal::order<DenseMIABase>::value> &_dims): MIA<DenseMIABase<Derived > >(_dims) {}

    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    bool operator==(const MIA<otherDerived>& otherMIA) const;

    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    bool fuzzy_equals(const MIA<otherDerived> & otherMIA,data_type precision ) const;



    //! Returns scalar data at given linear index
    inline const_data_type_ref atIdx(index_type idx) const{

        //return lin index
        return derived().atIdx(idx);
    }

    //! Returns scalar data at given linear index
    inline data_type_ref atIdx(index_type idx){

        //return lin index
        return derived().atIdx(idx);
    }








    data_iterator data_begin()
    {


        return derived().data_begin();
    }

    data_iterator data_end()
    {


        return derived().data_end();
    }

    const_data_iterator data_begin() const
    {


        return derived().data_begin();
    }

    const_data_iterator data_end() const
    {


        return derived().data_end();
    }




    template<class otherDerived>
    bool operator==(SparseMIABase<otherDerived>& otherMIA) const
    {
        return otherMIA==*this;
    }

    template<class otherDerived>
    bool operator!=(SparseMIABase<otherDerived>& otherMIA) const
    {
        return !(*this==otherMIA.derived());
    }

    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    bool operator!=(const MIA<otherDerived>& otherMIA) const
    {
        return !(*this==otherMIA.derived());
    }


    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type  plus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.plus_equal(b,index_order);
        return c;
    }


    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type  minus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.minus_equal(b,index_order);
        return c;
    }
    //!Component-wise addition with another MIA
    template<class otherDerived>
    typename MIAMergeReturnType<Derived,otherDerived>::type operator+(const MIA<otherDerived> &b)const{
        return this->outside_merge(b,std::plus<data_type>());
    }
    //!Component-wise subtraction with another MIA
    template<class otherDerived>
    typename MIAMergeReturnType<Derived,otherDerived>::type operator-(const MIA<otherDerived> &b)const{
        return this->outside_merge(b,std::minus<data_type>());
    }

    //! Add a single value to all elements (non-zero elements if sparse)
    typename MIANonlinearFuncType<Derived>::type operator+(data_type _data) const{
        typedef typename MIANonlinearFuncType<Derived>::type retType;
        retType ret(final_derived());

        for(auto it=ret.data_begin(); it< ret.data_end();++it)
            *it=*it+_data;
        return ret;
    }

    template<class otherDerived,typename index_param_type,typename Op>
    typename MIAMergeReturnType<Derived,otherDerived>::type  outside_merge(const MIA<otherDerived> &b,Op op,const std::array<index_param_type,mOrder>& index_order) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.merge(b,op,index_order);
        return c;
    }

    template<class otherDerived,typename Op>
    typename MIAMergeReturnType<Derived,otherDerived>::type  outside_merge(const MIA<otherDerived> &b,Op op) const{

        static_assert(mOrder==internal::order<otherDerived>::value,"Orders of two MIAs must be the same to perform addition");
        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.merge(b,op);
        return c;
    }





    template
    <
        size_t new_order=mOrder
    >
    const ImplicitMIA<data_type,new_order,true> view(std::array<Range<index_type>,mOrder> & ranges) const{
        const ImplicitMIA<data_type,new_order,true> ret=do_view<new_order>(ranges);
        return ret;
    }

    template
    <
        size_t new_order=mOrder
    >
    ImplicitMIA<data_type,new_order,true> view(std::array<Range<index_type>,mOrder> & ranges){
         return do_view<new_order>(ranges);
    }


    //! Returns a view (or subblock) of the MIA.
    /*!
        \param[in] ranges variadic parameter specifying range of the MIA view. Will assert a compile error if size of ranges!=mOrder or if any of Ranges datatype are not Range
        \return An ImplicitMIA that references *this's underlying raw data. Note changing the data in the ImplicitMIA will change *this. If you want a copy, call the make_explicit() function
                of the returned ImplicitMIA.
    */
    template<typename... Ranges>
    ImplicitMIA<data_type,internal::get_range_count<Ranges...>::range_count::value,true> view(Ranges... ranges) const {
        static_assert(internal::check_ranges<DenseMIABase,Ranges...>::type::value,"Range or integral datatypes must be passed to MIA when creating a view. Number of arguments must be equal to MIA order.");
        typedef typename internal::get_range_count<Ranges...>::range_count range_count;

        std::array<Range<index_type>,mOrder> range_array;
        internal::get_range_array(range_array.begin(),ranges...);
        return view<(size_t)(range_count::value)>(range_array);
    }






    //DenseMIABase<Derived>::contract_attract(const std::array<int,no_con_indices> & contract_indices,const std::array<int,no_con_partition> & contract_partition,const std::array<int,no_attract_indices> & attract_indices, const std::array<int,no_attract_partition> & attract_partition) const


//



//    template<class otherDerived,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
//    typename MIAMergeReturnType<Derived,otherDerived>::type  plus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{
//        std::plus<data_type> op;
//        return implicit_merge(b.derived(),op,index_order);
//    }
//
//    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
//    typename MIAMergeReturnType<Derived,otherDerived>::type  plus_(const MIA<otherDerived> &b) const{
//        std::plus<data_type> op;
//        return implicit_merge(b.derived(),op);
//    }
//
//    template<class otherDerived,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
//    typename MIAMergeReturnType<Derived,otherDerived>::type  minus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{
//        std::minus<data_type> op;
//        return implicit_merge(b.derived(),op,index_order);
//    }
//
//    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
//    typename MIAMergeReturnType<Derived,otherDerived>::type  minus_(const MIA<otherDerived> &b) const{
//        std::minus<data_type> op;
//        return implicit_merge(b.derived(),op);
//    }






    void print() const
    {
        print_array(this->dims(),"Dimensions");
        std::cout << "Index\t Data" << std::endl;
        for(size_t idx=0; idx<this->dimensionality();++idx){
            if(this->atIdx(idx))
                std::cout << idx << "\t " << this->atIdx(idx) << std::endl;
        }
        std::cout << std::endl;
    }


    //! Function for performing contraction and attraction. Best not to call directly, instead use the MIA algebra unary operation functionality, e.g., a(i,i)
    /*!
        \param[in] contract_indices list of indices undergoing a contraction
        \param[in] contract_partitions if more than one set of contractions is taking place, specifies how to partition contract_indices into corresponding sets of contractions
        \param[in] attract_indices list of indices undergoing an attraction
        \param[in] attract_partitions if more than one set of attractions is taking place, specifies how to partition attract_indices into corresponding sets of attractions
    */
    template<size_t no_con_indices,size_t no_con_partitions,size_t no_attract_indices,size_t no_attract_partitions>
    typename MIAUnaryType<Derived,internal::order<Derived>::value-no_con_indices-no_attract_indices+no_attract_partitions>::type contract_attract(const std::array<int,no_con_indices> & contract_indices,const std::array<int,no_con_partitions> & contract_partitions,const std::array<int,no_attract_indices> & attract_indices,const std::array<int,no_attract_partitions> & attract_partitions) const;




    //! Flattens the MIA to a Lattice by creating a copy of the data.
    /*!
        \param[in] row_indices indices to map to the lattice rows - will perserve ordering
        \param[in] column_indices indices to map to the lattice columns - will perserve ordering
        \param[in] tab_indices indices to map to the lattice tabs - will perserve ordering
        \return DenseLattice class that owns a copy of this's raw data
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    DenseLattice<data_type> toLatticeCopy(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const;


    DenseLattice<data_type> toStraightLatticeCopy(size_t number_of_row_indices,size_t number_of_column_indices) const;

protected:

     //! Routine for performing multiplication without the need for a lattice multiplication. Note, since b is const, multiplying by itself is fine
    template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
    noLatticeMult(const MIA<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx) const;

    //! Routine for performing multiplication without the need for a lattice multiplication
    template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
    noLatticeMult(SparseMIABase<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx) const{
            auto c= b.noLatticeMult(*this,r_inter_idx,r_outer_idx,l_inter_idx,l_outer_idx);
            auto c_reorder=internal::concat_index_arrays(internal::createAscendingIndex<L_outer>(R_outer),internal::createAscendingIndex<R_outer>(),internal::createAscendingIndex<R_inter>(L_outer+R_outer));
            //print_array(c.linIdxSequence(),"Old sort order");
            c.set_linIdxSequence(c_reorder);
            //print_array(c.linIdxSequence(),"New sort order");
            //print_array(c.dims(),"Old dims");
            auto new_c_dims=c.dims();
            internal::reorder_from(c.dims(), c_reorder,new_c_dims);
            //print_array(new_c_dims,"New dims");
            c.set_dims(new_c_dims);
            return c;

    }


    //! Returns scalar data at given indices
    /*!
        \param[in] ranges array of Range's specifying range of the MIA view. Will assert a compile error if size of ranges!=mOrder or if any of Ranges datatype are not Range
        \return An ImplicitMIA that references *this's underlying raw data. Note changing the data in the ImplicitMIA will change *this
    */
    template
    <
        size_t new_order=mOrder,
        typename boost::enable_if_c<
            new_order==internal::order<Derived>::value,
            int
        >::type=0
    >
    ImplicitMIA<data_type,new_order,true> do_view( std::array<Range<index_type>,mOrder> & ranges) const;

    template
    <
        size_t new_order=mOrder,
        typename boost::enable_if_c<
            new_order<internal::order<Derived>::value,
            int
        >::type=0
    >
    ImplicitMIA<data_type,new_order,true> do_view( std::array<Range<index_type>,mOrder> & ranges) const;






    //! Common routine for merge operations, such as add or subtract that result in an ImplicitMIA
    template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type
    implicit_merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order) const;


    //! Common routine for merge operations, such as add or subtract that result in an ImplicitMIA
    template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type
    implicit_merge(const MIA<otherDerived> &b,const Op& op) const;





private:

    template <class D1,class D2,bool D3,size_t D4> friend class MIA_Atom;
    template <class E1> friend class DenseMIABase;
    template <class F1> friend class MIA;


};




template<typename Derived>
template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
auto  DenseMIABase<Derived>::noLatticeMult(const MIA<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx) const
->typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
{




    return internal::implicitNoLatticeMult(this->final_derived(),b.final_derived(),l_inter_idx,l_outer_idx,r_inter_idx,r_outer_idx);


}

template<typename Derived>
template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
auto DenseMIABase<Derived>::toLatticeCopy(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const ->DenseLattice<data_type>
{

    return internal::latticeCopy(*this,row_indices,column_indices,tab_indices);




}



template<typename Derived>
auto DenseMIABase<Derived>::toStraightLatticeCopy(size_t number_of_row_indices, size_t number_of_column_indices) const ->DenseLattice<data_type>
{

    index_type row_size=std::accumulate(this->dims().begin(),this->dims().begin()+number_of_row_indices,1,std::multiplies<index_type>());
    index_type column_size=std::accumulate(this->dims().begin()+number_of_row_indices,this->dims().begin()+number_of_row_indices+number_of_column_indices,1,std::multiplies<index_type>());
    index_type tab_size=std::accumulate(this->dims().begin()+number_of_row_indices+number_of_column_indices,this->dims().end(),1,std::multiplies<index_type>());

    DenseLattice<data_type> lat(row_size,column_size,tab_size);
    for(size_t ctr=0;ctr<this->dimensionality();++ctr){
        lat.atIdx(ctr)=this->atIdx(ctr);
    }

    return lat;



}


template<typename Derived>
template<typename otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type >
bool DenseMIABase<Derived>::operator==(const MIA<otherDerived> & otherMIA ) const
{

    if(this->m_dims!=otherMIA.dims())
        return false;

    for(size_t idx=0;idx<this->dimensionality();++idx)
        if(this->atIdx(idx)!=otherMIA.atIdx(idx))
            return false;

    return true;

}


template<typename Derived>
template<typename otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
bool DenseMIABase<Derived>::fuzzy_equals(const MIA<otherDerived> & otherMIA,data_type precision ) const
{
    if(this->m_dims!=otherMIA.dims())
        return false;

    for(size_t idx=0;idx<this->dimensionality();++idx)
        if(!isEqualFuzzy(this->atIdx(idx),otherMIA.atIdx(idx),precision)){
            //std::cout << "Triggered " << *it1 << " " << *it2 << std::endl;
            return false;
        }

    return true;

}



template<typename Derived>
template
    <
        size_t new_order,
        typename boost::enable_if_c<
            new_order < internal::order<Derived>::value,
            int
        >::type
    >
auto DenseMIABase<Derived>::do_view(std::array<Range<index_type>,mOrder> &ranges)const->ImplicitMIA<data_type,new_order,true>  {
#ifdef LIBMIA_CHECK_DIMS
    size_t constant_count=0;
    for(auto &i:ranges){
        if(i.mEnd<i.mBegin)
            throw MIAParameterException("End of Range must be greater than beginning.");
        if((i.mEnd-i.mBegin)%i.mStep)
            throw MIAParameterException("Bounds of Range must be divisible by step size.");
        if(i.mEnd==i.mBegin)
            constant_count++;

    }
    if (constant_count>mOrder-new_order)
        throw MIAParameterException("Number of non-singleton ranges must be equal to the new_order template parameter.");
#endif

    typedef ImplicitMIA<data_type,new_order,true> retType; //the true tparam ensures that it returns references to the same raw data as *true
    //calculate returing dimensions
    std::array<index_type,new_order> retDims;
    //calculate returing dimensions
    std::array<size_t,new_order> retDimIndices;
    //holds the indices that remain constant (ie Range is only 1)
    std::array<index_type,mOrder> baseline_indices;
    size_t retIdx=0;
    for(size_t idx=0;idx<mOrder;++idx){
        if(ranges[idx].mEnd==-1)
            ranges[idx].mEnd=this->dim(idx);
        if(ranges[idx].mEnd>ranges[idx].mBegin){
            retDimIndices[retIdx]=idx;
            retDims[retIdx++]=(ranges[idx].mEnd-ranges[idx].mBegin)/ranges[idx].mStep;
        }
        else if(ranges[idx].mEnd==ranges[idx].mBegin)
            baseline_indices[idx]=ranges[idx].mBegin; //store the constant-valued index
    }
    //now create a lambda function that returns a reference to this's data based on the ranges provided
    typedef typename internal::function_type<retType>::type function_type;

    function_type func=[this,ranges,retDims,retDimIndices,baseline_indices](index_type idx)->data_type&{
        //std::cout << "idx " << idx <<std::endl;

        auto indices=internal::ind2sub(idx,retDims); //get the indices in the subview
        //print_array(indices,"indices");
        auto expanded_indices=baseline_indices;
        //convert the indices to the indices in *this
        for(size_t idx=0;idx<new_order;++idx){
            expanded_indices[retDimIndices[idx]]=indices[idx]*ranges[retDimIndices[idx]].mStep+ranges[retDimIndices[idx]].mBegin;
        }

        //print_array(indices,"indices");
        //the const cast is ugly, but we're relying on const overloading for the view function and also for the ImplicitMIA to ensure any const references are not modified
        //a less hacky solution would be nicer
        return const_cast<data_type&>(this->at(expanded_indices)); //return the corresponding data from *this
    };

    return retType(func,retDims);
}


template<typename Derived>
template
    <
        size_t new_order,
        typename boost::enable_if_c<
            new_order==internal::order<Derived>::value,
            int
        >::type
    >
auto DenseMIABase<Derived>::do_view(std::array<Range<index_type>,mOrder> &ranges) const->ImplicitMIA<data_type,new_order,true> {
#ifdef LIBMIA_CHECK_DIMS
    for(auto &i:ranges){
        if(i.mEnd<i.mBegin)
            throw MIAParameterException("End of Range must be greater than beginning.");
        if((i.mEnd-i.mBegin)%i.mStep)
            throw MIAParameterException("Bounds of Range must be divisible by step size.");

    }
#endif

    typedef ImplicitMIA<data_type,new_order,true> retType; //the true tparam ensures that it returns references to the same raw data as *true
    //calculate returing dimensions
    std::array<index_type,new_order> retDims;
    size_t retIdx=0;
    for(size_t idx=0;idx<new_order;++idx){
        if(ranges[idx].mEnd==-1)
            ranges[idx].mEnd=this->dim(idx);
        retDims[retIdx++]=(ranges[idx].mEnd-ranges[idx].mBegin)/ranges[idx].mStep;
    }
    //now create a lambda function that returns a reference to this's data based on the ranges provided
    typedef typename internal::function_type<retType>::type function_type;

    function_type func=[this,ranges,retDims](index_type idx)->data_type&{
        //std::cout << "idx " << idx <<std::endl;

        auto indices=internal::ind2sub(idx,retDims); //get the indices in the subview
        //print_array(indices,"indices");
        //convert the indices to the indices in *this
        for(size_t idx=0;idx<new_order;++idx){
            indices[idx]=indices[idx]*ranges[idx].mStep+ranges[idx].mBegin;
        }

        //the const cast is ugly, but we're relying on const overloading for the view function and also for the ImplicitMIA to ensure any const references are not modified
        //a less hacky solution would be nicer
        return const_cast<data_type&>(this->at(indices)); //return the corresponding data from *this
    };

    return retType(func,retDims);
}

template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
typename MIAMergeReturnType<Derived,otherDerived>::type
DenseMIABase<Derived>::implicit_merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order) const
{


    //std::cout << "Implicit" << std::endl;
    #ifdef LIBMIA_CHECK_DIMS
    this->check_merge_dims(b,index_order);
    #endif
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");


    return internal::perform_implicit_merge(*this, b,op,index_order);


}

template<typename Derived>
template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
typename MIAMergeReturnType<Derived,otherDerived>::type
DenseMIABase<Derived>::implicit_merge(const MIA<otherDerived> &b,const Op& op) const
{



    #ifdef LIBMIA_CHECK_DIMS
    this->check_merge_dims(b);
    #endif



    return internal::perform_implicit_merge(*this, b,op);


}


template<typename Derived>
template<size_t no_con_indices,size_t no_con_partition,size_t no_attract_indices, size_t no_attract_partition>
typename MIAUnaryType<Derived,internal::order<Derived>::value-no_con_indices-no_attract_indices+no_attract_partition>::type
DenseMIABase<Derived>::contract_attract(const std::array<int,no_con_indices> & contract_indices,const std::array<int,no_con_partition> & contract_partition,const std::array<int,no_attract_indices> & attract_indices, const std::array<int,no_attract_partition> & attract_partition) const
{

    static_assert(no_con_indices+no_attract_indices<=mOrder,"Number of indices specified for contraction and/or attraction must not exceed mOrder");
    #ifdef LIBMIA_CHECK_DIMS
    //this->check_contract_indices(contract_indices,attract_indices);
    #endif
    typedef typename MIAUnaryType<Derived,internal::order<Derived>::value-no_con_indices-no_attract_indices+no_attract_partition>::type retType;

    //extract indices not undergoing a contraction or attraction
    auto copy_contract=internal::concat_index_arrays(contract_indices,attract_indices);
    std::sort(copy_contract.begin(),copy_contract.end());
    auto other_indices=internal::get_remaining_indices<size_t,no_con_indices+no_attract_indices,mOrder>(copy_contract);
    //get their dimensionality
    std::array<index_type,other_indices.size()> otherDims;
    auto other_dimensionality=internal::reorder_from(this->dims(), other_indices,otherDims);


    //print_array(contract_indices,"contract_indices");
    //print_array(contract_partition,"contract_partition");

//    print_array(attract_indices,"attract_indices");
//    print_array(attract_partition,"attract_partition");

    //print_array(other_indices,"other_indices");
    //

    //get size of contraction indices in each set of contractions (each range within the same partition or set should be all identical)
    std::array<size_t,no_con_partition> contract_index_ranges;
    size_t cur_idx=0;
    for(size_t i=0;i<no_con_partition;++i)
    {
       contract_index_ranges[i]= this->dim(contract_indices[cur_idx]);
       cur_idx+=contract_partition[i];

    }
    //get size of attracion index ranges in each set of attractions (each range within the same partition or set should be all identical)
    std::array<size_t,no_attract_partition> attract_index_ranges;
    cur_idx=0;
    size_t attract_dimensionality=1;
    for(size_t i=0;i<no_attract_partition;++i)
    {
       attract_index_ranges[i]= this->dim(attract_indices[cur_idx]);
       attract_dimensionality*=attract_index_ranges[i];
       cur_idx+=attract_partition[i];

    }
    //add the attraction index ranges to the otherDims to get the returning dimensionality
    auto retDims=internal::concat_index_arrays(otherDims,attract_index_ranges);

    //print_array(other_indices,"other_indices");
    //print_array(retDims,"retDims");

    retType ret(retDims);

    //loop through all index locations of the returning MIA not undergoing an attraction or contraction

    for(size_t i=0;i<other_dimensionality;++i){
        //calculate the current index location in the original MIA
        auto other_i_idx=internal::sub2ind(internal::ind2sub(i,otherDims), other_indices, this->dims()); //get location of current index in the original MIA

        //for each attraction index location (if none, attract_dimensionality will be 1), calculate the corresponding element value
        for(size_t j=0;j<attract_dimensionality;++j)
        {

            size_t other_j_idx=0;
            //get current attract index values in returning MIA (could be more than one index value if there is more than one set of attractions)
            auto attract_idx=internal::ind2sub(j,attract_index_ranges);
            auto cur_partition_begin=attract_indices.begin();
            for(size_t k=0;k<no_attract_partition;++k){
                //get the index location of current attraction index in original MIA. E.g., if n, will be {n,n} or {n,n,} for some set of indices in original MIA
                other_j_idx+=internal::get_contract_idx(attract_idx[k], cur_partition_begin,cur_partition_begin+attract_partition[k], this->dims());
                cur_partition_begin+=attract_partition[k];
            }
            //caculate returing MIAs data value. If a contraction is taking place, a contraction will be performed. Otherwise, it just returns the appropriate data value in the original MIA
            ret.atIdx(i+j*other_dimensionality)=internal::collect_contract_partitions<no_con_partition,FinalDerived,decltype(contract_indices.end()),no_con_partition>(this->final_derived(),other_i_idx+other_j_idx,contract_index_ranges,contract_indices.end(),contract_partition);
        }


    }


    return ret;

}



//template<typename Derived>
//template<class otherDerived,typename index_param_type>
//auto DenseMIABase<Derived>::plus_equal(const DenseMIABase<otherDerived> &b,const std::array<index_param_type,DenseMIABase::order>& index_order)->DenseMIABase &
//{
//
//    std::plus<data_type> op;
//    merge(b,op,index_order);
//    return *this;
//
//}
//
//template<typename Derived>
//template<class otherDerived,typename index_param_type>
//auto DenseMIABase<Derived>::minus_equal(const DenseMIABase<otherDerived> &b,const std::array<index_param_type,DenseMIABase::order>& index_order)->DenseMIABase &
//{
//
//    std::minus<data_type> op;
//    merge(b,op,index_order);
//    return *this;
//
//}


/*! @} */

}






#endif // DENSEMIABASE_H_INCLUDED
