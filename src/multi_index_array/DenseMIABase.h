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
    inline const_data_type_ref atIdx(index_type idx) const{

        //return lin index
        return derived().atIdx(idx);
    }

    //! Returns scalar data at given linear index
    inline data_type_ref atIdx(index_type idx){

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


    DenseLattice<data_type> toStraightLatticeCopy(size_t number_of_row_indices,size_t number_of_column_indices) const;



    //! Common routine for merge operations, such as add or subtract that result in an ImplicitMIA
    template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type
    implicit_merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order) const;


    //! Common routine for merge operations, such as add or subtract that result in an ImplicitMIA
    template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    typename MIAMergeReturnType<Derived,otherDerived>::type
    implicit_merge(const MIA<otherDerived> &b,const Op& op) const;

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


    template<class otherDerived,typename index_param_type,typename Op>
    typename MIAMergeReturnType<Derived,otherDerived>::type  outside_merge(const MIA<otherDerived> &b,Op op,const std::array<index_param_type,mOrder>& index_order) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.merge(b,op,index_order);
        return c;
    }

    template<class otherDerived,typename Op>
    typename MIAMergeReturnType<Derived,otherDerived>::type  outside_merge(const MIA<otherDerived> &b,Op op) const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        CType c(*this);
        c.merge(b,op);
        return c;
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





protected:

    SolveInfo mSolveInfo;


private:



};




template<typename Derived>
template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
auto  DenseMIABase<Derived>::noLatticeMult(const MIA<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
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
    internal::reorder_from(this->dims(), other_indices,otherDims);


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
    size_t other_dimensionality=ret.dimensionality()/attract_dimensionality;
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
