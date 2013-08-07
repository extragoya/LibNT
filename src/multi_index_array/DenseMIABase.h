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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/dynamic_bitset.hpp>

#include "Util.h"
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
    DenseMIABase(Dims... dims): MIA<DenseMIABase<Derived > >(dims...),mSolveInfo(NoInfo) {}


    DenseMIABase(): MIA<DenseMIABase<Derived > >(),mSolveInfo(NoInfo) {}


    DenseMIABase(std::array<index_type,internal::order<DenseMIABase>::value> &_dims): MIA<DenseMIABase<Derived > >(_dims),mSolveInfo(NoInfo) {}

    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    bool operator==(const MIA<otherDerived>& otherMIA) const;

    template<class otherDerived,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    bool fuzzy_equals(const MIA<otherDerived> & otherMIA,data_type precision ) const;



    //! Returns scalar data at given linear index
    const_data_type_ref atIdx(index_type idx) const{

        //return lin index
        return derived().atIdx(idx);
    }

    //! Returns scalar data at given linear index
    data_type_ref atIdx(index_type idx){

        //return lin index
        return derived().atIdx(idx);
    }

    SolveInfo solveInfo() const{
        return mSolveInfo;
    }

    void setSolveInfo(SolveInfo _solveInfo){
        mSolveInfo=_solveInfo;
    }



    //! Flattens the MIA to a Lattice by creating a copy of the data.
    /*!
        \param[in] row_indices indices to map to the lattice rows - will perserve ordering
        \param[in] column_indices indices to map to the lattice columns - will perserve ordering
        \param[in] tab_indices indices to map to the lattice tabs - will perserve ordering
        \return DenseLattice class that owns a copy of this's raw data
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    DenseLattice<data_type> toLatticeCopy(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const;






    //! Common routine for merge operations, such as add or subtract that result in an ImplicitMIA
    template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type
    implicit_merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order) const;

    //! Routine for performing multiplication without the need for a lattice multiplication
    template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
    noLatticeMult(const MIA<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

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
    typename MIAMergeReturnType<Derived,otherDerived>::type  plus_(const SparseMIABase<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.plus_equal(b,index_order);
        return c;
    }

    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type  minus_(const SparseMIABase<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.minus_equal(b,index_order);
        return c;
    }

    template<class otherDerived,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type  plus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{
        std::plus<data_type> op;
        return implicit_merge(b.derived(),op,index_order);
    }

    template<class otherDerived,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type  minus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{
        std::minus<data_type> op;
        return implicit_merge(b.derived(),op,index_order);
    }





    void print() const
    {
        std::cout << "Index\t Data" << std::endl;
        for(size_t idx=0; idx<this->dimensionality();++idx){
            if(this->atIdx(idx))
                std::cout << idx << "\t " << this->atIdx(idx) << std::endl;
        }
        std::cout << std::endl;
    }


    template<typename index_param_type>
    void inplace_permute(const std::array<index_param_type,internal::order<Derived>::value> & reshuffle_order);


    //!
    /*!
        Based on An Optimal Index Reshuffle Algorithm for Multidimensional Arrays and Its Applications for Parallel Architectures by Chris H.Q. Ding and the modification
        in Sec. III A by Jie et al.'s article: A High Efficient In-place Transposition Scheme for Multidimensional Arrays
    */
    template<typename... Dims>
    void inplace_permute(Dims... dims){
        static_assert(internal::check_mia_constructor<DenseMIABase,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");

        std::array<index_type, mOrder> reshuffle_order{{dims...}};
        inplace_permute(reshuffle_order);

    }



protected:

    SolveInfo mSolveInfo;


private:



};

template<typename Derived>
template<typename index_param_type>
void DenseMIABase<Derived>::inplace_permute(const std::array<index_param_type,internal::order<Derived>::value>& reshuffle_order){

    //first check that the reshuffle_order array isn't just {0,1,...mOrder}
    //if it is, we do nothing
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

    size_t check;
    for(check=0;check<reshuffle_order.size();++check){
        if(reshuffle_order[check]!=(index_param_type)check)
            break;
    }
    if (check==reshuffle_order.size())
        return;

    boost::dynamic_bitset<> bit_array(this->dimensionality()); // all 0's by default
    std::array<index_type,mOrder> new_dims; //stores new dimensions
    internal::reorder_from(this->dims(),reshuffle_order,new_dims); //get new dims
    index_type ioffset;

    std::array<index_type,mOrder> dim_accumulator; //precompute the demoninators needed to convert from linIdx to a full index, using new_dims
    for(size_t i=0;i<mOrder;++i){
        dim_accumulator[i]=std::accumulate(new_dims.begin(),new_dims.begin()+i,1,std::multiplies<index_type>());
    }

    dim_accumulator=internal::reorder_to(dim_accumulator,reshuffle_order); //reorder the denominators based on the reshuffle order
    //create a function that converts from a linIdx of the permuted array to a linIdx of the old array (see Jie et al.'s article: A High Efficient In-place Transposition Scheme for Multidimensional Arrays)
    //EDIT - using the lambda slowed down the implementation, so it's just hard-coded now
//    auto func=[this,&dim_accumulator](const index_type from_lin_idx){
//        index_type to_lin_idx=0;
//        index_type multiplier=1;
//        for(size_t i=0;i<mOrder;++i){
//            to_lin_idx+=(from_lin_idx/dim_accumulator[i])%this->dim(i)*multiplier; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly
//            multiplier*=this->dim(i);
//        }
//        return to_lin_idx;
//    };
    index_type touch_ctr=0;
    //iterate through the entire bit array
    for(index_type start_idx=0;start_idx<this->dimensionality();++start_idx){
        //if we've found a location that hasn't been touched, we've found a new vacancy cycle
        if(bit_array[start_idx]==0){

            auto temp=this->atIdx(start_idx); //get the start of the cycle
            ioffset=start_idx; //permuted array linIdx
            while(true){

                index_type ioffset_next=0;
                index_type multiplier=1;
                for(size_t i=0;i<mOrder;++i){
                    ioffset_next+=(ioffset/dim_accumulator[i])%this->dim(i)*multiplier; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly
                    multiplier*=this->dim(i);
                }
                //ioffset_next=func(ioffset); //get the location of ioffset in the old array


                bit_array[ioffset]=1; //touch the current data location
                ++touch_ctr;
                if(ioffset_next==start_idx){ //if we've cycled to the start, then finish the cycle and break
                    if(ioffset_next!=ioffset)
                        this->atIdx(ioffset)=temp;
                    break;
                }

                this->atIdx(ioffset)=this->atIdx(ioffset_next); //set the data element in the permuted array to its location in the old array
                ioffset=ioffset_next; //update the location in the current cycle

            }
            if (touch_ctr==this->dimensionality()){
                break;
            }

        }


    }
    this->m_dims=new_dims;


}


template<typename Derived>
template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
auto  DenseMIABase<Derived>::noLatticeMult(const MIA<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
{




    return internal::implicitNoLatticeMult(this->final_derived(),b.final_derived(),l_inter_idx,l_outer_idx,r_inter_idx,r_outer_idx);
//    //statically check number of indices match up
//    static_assert(internal::check_index_compatibility<index_type,array_type>::type::value,"Must use an array convertable to index_type");
//    typedef typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type RetType; //get return type
//    typedef typename internal::index_type<otherDerived>::type b_index_type;
//
//    std::array<index_type, L_inter> l_inter_dims;
//    std::array<index_type, L_outer> l_outer_dims;
//    std::array<b_index_type, R_inter> r_inter_dims;
//    std::array<b_index_type, R_outer> r_outer_dims;
//    std::array<typename internal::index_type<RetType>::type,L_outer+R_outer+L_inter> return_dims;
//
//
//    //get inter and outer dimensionality and the individual dimensions that make up that number - should default to one if any of the arrays are empty
//    size_t l_inter_size=internal::reorder_from(this->dims(), l_inter_idx,l_inter_dims);
//    internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
//    if(l_inter_dims.size()!=r_inter_dims.size() || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
//        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");
//    size_t l_outer_size=internal::reorder_from(this->dims(), l_outer_idx,l_outer_dims);
//    size_t r_outer_size=internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);
//
//
//    //concatenate the pulled dimensions into a consolidated set of dimensions for the return mia
//
//    auto temp_end=std::copy(l_outer_dims.begin(),l_outer_dims.end(),return_dims.begin());
//    temp_end=std::copy(r_outer_dims.begin(),r_outer_dims.end(),temp_end);
//    temp_end=std::copy(l_inter_dims.begin(),l_inter_dims.end(),temp_end);
//    RetType C(return_dims);
//
//
//    index_type l_inter_linear_idx=0, l_outer_linear_idx=0;
//    b_index_type r_outer_linear_idx=0,r_inter_linear_idx=0;
//    for (index_type k=0;k<l_inter_size;++k){ //loop through every set of element-wise indices (in linearized form)
//        //compute the offset of the current elementwise index for both operands
//        l_inter_linear_idx=internal::sub2ind(internal::ind2sub(k,l_inter_dims),l_inter_idx,this->dims());
//        r_inter_linear_idx=internal::sub2ind(internal::ind2sub(k,r_inter_dims),r_inter_idx,b.dims());
//        //now loop through every possible set of outer-wise indices in both operands (these are not shared)
//        for (index_type j=0;j<l_outer_size;++j){
//            l_outer_linear_idx=l_inter_linear_idx+internal::sub2ind(internal::ind2sub(j,l_outer_dims),l_outer_idx,this->dims()); //compute offset
//            for (index_type i=0;i<r_outer_size;++i){
//                r_outer_linear_idx=r_inter_linear_idx+internal::sub2ind(internal::ind2sub(i,r_outer_dims),r_outer_idx,b.dims()); //compute offset
//                //set the element at the appropriate offset in C to the product of the two operand elements
//                C.atIdx(j+i*l_outer_size+k*r_outer_size*l_outer_size)=this->atIdx(l_outer_linear_idx)*b.atIdx(r_outer_linear_idx);
//
//            }
//        }
//
//
//    }
//    return C;

}

template<typename Derived>
template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
auto DenseMIABase<Derived>::toLatticeCopy(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const ->DenseLattice<data_type>
{

    return internal::latticeCopy(*this,row_indices,column_indices,tab_indices);

//    static_assert(internal::check_index_compatibility<index_type,idx_typeR>::type::value,"Must use an array convertable to index_type");
//    static_assert(internal::check_index_compatibility<index_type,idx_typeC>::type::value,"Must use an array convertable to index_type");
//    static_assert(internal::check_index_compatibility<index_type,idx_typeT>::type::value,"Must use an array convertable to index_type");
//    static_assert(R+C+T==mOrder,"Size of all three arrays must equal mOrder");
//    //statically check number of indices match up
//    size_t row_size=1, column_size=1, tab_size=1;
//    std::array<index_type, R> row_dims;
//    std::array<index_type, C> column_dims;
//    std::array<index_type, T> tab_dims;
//
//    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
//    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";
//
//    row_size=internal::reorder_from(this->dims(), row_indices,row_dims);
//    column_size=internal::reorder_from(this->dims(), column_indices,column_dims);
//    tab_size=internal::reorder_from(this->dims(), tab_indices,tab_dims);
//
//    //std::cout<< "Tab dims " << tab_dims[0] << "\n";
//    DenseLattice<data_type> lat(row_size, column_size, tab_size);
//
//
//
//    index_type row_idx=0, column_idx=0,tab_idx=0;
//    for (index_type k=0;k<tab_size;++k){
//        tab_idx=internal::sub2ind(internal::ind2sub(k,tab_dims),tab_indices,this->m_dims);
//        for (index_type j=0;j<column_size;++j){
//            column_idx=tab_idx+internal::sub2ind(internal::ind2sub(j,column_dims),column_indices,this->m_dims);
//            for (index_type i=0;i<row_size;++i){
//                row_idx=column_idx+internal::sub2ind(internal::ind2sub(i,row_dims),row_indices,this->m_dims);
//
//                lat(i,j,k)=this->atIdx(row_idx);
//
//            }
//        }
//
//
//    }
//
//    //lat.print();
//    return lat;


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
template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
typename MIAMergeReturnType<Derived,otherDerived>::type
DenseMIABase<Derived>::implicit_merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order) const
{


    //std::cout << "Implicit" << std::endl;
    this->check_merge_dims(b,index_order);
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

//    typedef typename MIAMergeReturnType<Derived,otherDerived>::type retMIAType;
//    typedef typename internal::function_type<retMIAType>::type function_type;
//    function_type _function;
//
//    if(internal::check_ascending(index_order)){ //if there is no index shuffle, the function is simpler and faster
//
//       _function=[this,&b,op](index_type idx){
//            return op(this->atIdx(idx),this->derived().convert(b.atIdx(idx)));
//            //return this->atIdx(idx)+b.atIdx(idx);
//        };
//    }
//    else{
//        _function=[this,&b,op,index_order](index_type idx){
//            return op(this->atIdx(idx),this->derived().convert(b.atIdx(internal::sub2ind(this->ind2sub(idx),index_order,b.dims()))));
//        };
//    }
////    function_type _function=[](index_type idx){
////        return 0;
////        //return op(this->atIdx(idx),b.atIdx(internal::sub2ind(this->ind2sub(idx),index_order,b.dims())));
////    };
//
//    return retMIAType(_function,this->dims());
    return internal::perform_implicit_merge(*this, b,op,index_order);


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
