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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/mpl/apply.hpp>

#include "Util.h"
#include "MIA.h"
#include "DenseLattice.h"



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
struct index_type<DenseMIABase<Derived> >: public index_type<Derived> {};

template<class Derived>
struct order<DenseMIABase<Derived> >: public order<Derived> {};

template<class Derived>
struct data_iterator<DenseMIABase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct const_data_iterator<DenseMIABase<Derived> >: public const_data_iterator<Derived> {};

template<class Derived>
struct storage_iterator<DenseMIABase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<DenseMIABase<Derived> >: public const_storage_iterator<Derived> {};

}



//!  Base class for dense multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  dense multi-index array types. Provides operations and functions common to all dense
  multi-index arrays.

  \tparam Derived   should only be a dense multi-index array class.
*/
template <class Derived>
class DenseMIABase: public MIA<DenseMIABase<Derived > >
{
public:


    typedef typename internal::data_type<Derived>::type data_type;
    typedef typename internal::index_type<Derived>::type index_type;
    typedef typename internal::Data<Derived>::type Data;
    typedef typename internal::data_iterator<Derived>::type data_iterator;
    typedef typename internal::const_data_iterator<Derived>::type const_data_iterator;
    typedef typename internal::storage_iterator<Derived>::type storage_iterator;
    typedef typename internal::const_storage_iterator<Derived>::type const_storage_iterator;
    constexpr static size_t mOrder=internal::order<Derived>::value;
    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }
    /** \returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

    template<typename... Dims>
    DenseMIABase(Dims... dims): MIA<DenseMIABase<Derived > >(dims...),mSolveInfo(NoInfo) {}


    DenseMIABase(): MIA<DenseMIABase<Derived > >(),mSolveInfo(NoInfo) {}


    DenseMIABase(std::array<index_type,internal::order<DenseMIABase>::value> &_dims): MIA<DenseMIABase<Derived > >(_dims),mSolveInfo(NoInfo) {}

    template<class otherDerived>
    bool operator==(const DenseMIABase<otherDerived>& otherMIA);

    template<class otherDerived>
    bool fuzzy_equals(const DenseMIABase<otherDerived> & otherMIA,data_type precision );



    //!Assignment operator. Will call Derived's operator
    template<class otherDerived>
    DenseMIABase& operator=(const DenseMIABase<otherDerived>& otherMIA){
        return derived()=otherMIA;
    }


    //!Assignment operator. Will call Derived's operator
    DenseMIABase& operator=(const DenseMIABase& otherMIA){
        return derived()=otherMIA.derived();
    }

    //!  Assignment based on given order. Will call Derived's operator
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up

        \param[in] otherMIA the other MIA
        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {3,1,2} this->at(1,2,3)==otherMIA.at(2,3,1).
                                Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,class index_param_type>
    void assign(const MIA<otherDerived>& otherMIA,const std::array<index_param_type,internal::order<DenseMIABase>::value>& index_order){
        derived().assign(otherMIA.derived(),index_order);
    }



    //! Returns scalar data at given linear index
    const data_type& atIdx(index_type idx) const{

        //return lin index
        return *(derived().data_begin()+idx);
    }

    //! Returns scalar data at given linear index
    data_type& atIdx(index_type idx) {

        //return lin index
        return *(derived().data_begin()+idx);
    }

    SolveInfo solveInfo() const{
        return mSolveInfo;
    }

    void setSolveInfo(SolveInfo _solveInfo){
        mSolveInfo=_solveInfo;
    }

    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom.
    /*!
        For DenseMIAs, this function calls toLatticeCopy
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    DenseLattice<data_type> toLatticeExpression(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const
    {
        return toLatticeCopy(row_indices, column_indices, tab_indices);

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


    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
    template<typename otherDerived, typename Op,typename index_param_type>
    void  merge(const DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order);

    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
    template<typename otherDerived, typename Op,typename index_param_type>
    void  merge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order);

    //! Routine for performing multiplication without the need for a lattice multiplication
    template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
    noLatticeMult(const DenseMIABase<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

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


    storage_iterator begin()
    {


        return derived().begin();
    }

    storage_iterator end()
    {


        return derived().end();
    }

    const_storage_iterator begin() const
    {


        return derived().begin();
    }

    const_storage_iterator end() const
    {


        return derived().end();
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

    template<class otherDerived>
    bool operator!=(const DenseMIABase<otherDerived>& otherMIA) const
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

    //! Performs destructive add (+=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is added like: this->at(x,y,z)+=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,typename index_param_type>
    DenseMIABase & plus_equal(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){

        std::plus<data_type> op;
        merge(b.derived(),op,index_order);
        return *this;
    }


    //! Performs destructive subtract (-=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is subtracted like: this->at(x,y,z)-=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,typename index_param_type>
    DenseMIABase & minus_equal(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){
        std::minus<data_type> op;
        merge(b.derived(),op,index_order);
        return *this;
    }

    void print()
    {
        std::cout << "Index\t Data" << std::endl;
        for(auto it=this->data_begin();it<this->data_end();++it){
            if(*it)
                std::cout << it-this->data_begin() << "\t " << *it << std::endl;
        }
        std::cout << std::endl;
    }

protected:

    SolveInfo mSolveInfo;


private:



};

template<typename Derived>
template<typename otherDerived, typename array_type,size_t L_inter,size_t L_outer,size_t R_inter,size_t R_outer>
auto  DenseMIABase<Derived>::noLatticeMult(const DenseMIABase<otherDerived> &b,const std::array<array_type,L_inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,R_inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type
{
    //statically check number of indices match up
    static_assert(internal::check_index_compatibility<index_type,array_type>::type::value,"Must use an array convertable to index_type");
    typedef typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+L_inter>::type RetType; //get return type
    typedef typename internal::index_type<otherDerived>::type b_index_type;

    std::array<index_type, L_inter> l_inter_dims;
    std::array<index_type, L_outer> l_outer_dims;
    std::array<b_index_type, R_inter> r_inter_dims;
    std::array<b_index_type, R_outer> r_outer_dims;
    std::array<typename internal::index_type<RetType>::type,L_outer+R_outer+L_inter> return_dims;


    //get inter and outer dimensionality and the individual dimensions that make up that number - should default to one if any of the arrays are empty
    size_t l_inter_size=internal::reorder_from(this->dims(), l_inter_idx,l_inter_dims);
    internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
    if(l_inter_dims.size()!=r_inter_dims.size() || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");
    size_t l_outer_size=internal::reorder_from(this->dims(), l_outer_idx,l_outer_dims);
    size_t r_outer_size=internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);


    //concatenate the pulled dimensions into a consolidated set of dimensions for the return mia

    auto temp_end=std::copy(l_outer_dims.begin(),l_outer_dims.end(),return_dims.begin());
    temp_end=std::copy(r_outer_dims.begin(),r_outer_dims.end(),temp_end);
    temp_end=std::copy(l_inter_dims.begin(),l_inter_dims.end(),temp_end);
    RetType C(return_dims);


    index_type l_inter_linear_idx=0, l_outer_linear_idx=0;
    b_index_type r_outer_linear_idx=0,r_inter_linear_idx=0;
    for (index_type k=0;k<l_inter_size;++k){ //loop through every set of element-wise indices (in linearized form)
        //compute the offset of the current elementwise index for both operands
        l_inter_linear_idx=internal::sub2ind(internal::ind2sub(k,l_inter_dims),l_inter_idx,this->dims());
        r_inter_linear_idx=internal::sub2ind(internal::ind2sub(k,r_inter_dims),r_inter_idx,b.dims());
        //now loop through every possible set of outer-wise indices in both operands (these are not shared)
        for (index_type j=0;j<l_outer_size;++j){
            l_outer_linear_idx=l_inter_linear_idx+internal::sub2ind(internal::ind2sub(j,l_outer_dims),l_outer_idx,this->dims()); //compute offset
            for (index_type i=0;i<r_outer_size;++i){
                r_outer_linear_idx=r_inter_linear_idx+internal::sub2ind(internal::ind2sub(i,r_outer_dims),r_outer_idx,b.dims()); //compute offset
                //set the element at the appropriate offset in C to the product of the two operand elements
                C.atIdx(j+i*l_outer_size+k*r_outer_size*l_outer_size)=this->atIdx(l_outer_linear_idx)*b.atIdx(r_outer_linear_idx);

            }
        }


    }
    return C;

}

template<typename Derived>
template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
auto DenseMIABase<Derived>::toLatticeCopy(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const ->DenseLattice<data_type>
{

    static_assert(internal::check_index_compatibility<index_type,idx_typeR>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeC>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeT>::type::value,"Must use an array convertable to index_type");

    //statically check number of indices match up
    size_t row_size=1, column_size=1, tab_size=1;
    std::array<index_type, R> row_dims;
    std::array<index_type, C> column_dims;
    std::array<index_type, T> tab_dims;

    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    row_size=internal::reorder_from(this->dims(), row_indices,row_dims);
    column_size=internal::reorder_from(this->dims(), column_indices,column_dims);
    tab_size=internal::reorder_from(this->dims(), tab_indices,tab_dims);

    //std::cout<< "Tab dims " << tab_dims[0] << "\n";
    DenseLattice<data_type> lat(row_size, column_size, tab_size);



    index_type row_idx=0, column_idx=0,tab_idx=0;
    for (index_type k=0;k<tab_size;++k){
        tab_idx=internal::sub2ind(internal::ind2sub(k,tab_dims),tab_indices,this->m_dims);
        for (index_type j=0;j<column_size;++j){
            column_idx=tab_idx+internal::sub2ind(internal::ind2sub(j,column_dims),column_indices,this->m_dims);
            for (index_type i=0;i<row_size;++i){
                row_idx=column_idx+internal::sub2ind(internal::ind2sub(i,row_dims),row_indices,this->m_dims);

                lat(i,j,k)=this->atIdx(row_idx);

            }
        }


    }

    //lat.print();
    return lat;


}

template<typename Derived>
template<typename otherDerived>
bool DenseMIABase<Derived>::operator==(const DenseMIABase<otherDerived> & otherMIA )
{
    if(this->m_dims!=otherMIA.dims())
        return false;

    for(auto it1=this->data_begin(),it2=otherMIA.data_begin();it1<this->data_end();++it1,++it2)
        if(*it1!=*it2)
            return false;

    return true;

}


template<typename Derived>
template<typename otherDerived>
bool DenseMIABase<Derived>::fuzzy_equals(const DenseMIABase<otherDerived> & otherMIA,data_type precision )
{
    if(this->m_dims!=otherMIA.dims())
        return false;

    for(auto it1=this->data_begin(),it2=otherMIA.data_begin();it1<this->data_end();++it1,++it2)
        if(!isEqualFuzzy(*it1,*it2,precision)){
            std::cout << "Triggered " << *it1 << " " << *it2 << std::endl;
            return false;
        }

    return true;

}

template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type>
void  DenseMIABase<Derived>::merge(const DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order)
{

    this->check_merge_dims(b,index_order);
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");





    for(index_type curIdx=0; curIdx<this->m_dimensionality; ++curIdx)
    {
        this->atIdx(curIdx)=op(this->atIdx(curIdx),derived().convert(b.atIdx(internal::sub2ind(this->ind2sub(curIdx),index_order,b.dims()))));

    }


}

template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type>
void  DenseMIABase<Derived>::merge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::mOrder>& index_order)
{

    this->check_merge_dims(b,index_order);
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

    for(auto it=b.storage_begin();it<b.storage_end();++it){
        //the index values of b may correspond to a shuffled version of the default linear index
        auto default_order_idx=b.convert_to_default_sort(b.index_val(*it));
        //calculate the lhs_index based on how the MIA indices matched up
        auto lhs_index=internal::sub2ind_reorder(b.ind2sub(default_order_idx),index_order,b.dims());
        this->atIdx(lhs_index)=op(this->atIdx(lhs_index),derived().convert(b.data_val(*it)));

    }

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
