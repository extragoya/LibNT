// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef SPARSEMIABASE_H_INCLUDED
#define SPARSEMIABASE_H_INCLUDED

#include "LibMiaException.h"
#include "Util.h"
#include "IndexUtil.h"
#include "MIA.h"



//\defgroup
namespace LibMIA
{

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */

namespace internal
{



template<class Derived>
struct data_type<SparseMIABase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<SparseMIABase<Derived> >: public index_type<Derived> {};

template<class Derived>
struct order<SparseMIABase<Derived> >: public order<Derived> {};

template<class Derived>
struct data_iterator<SparseMIABase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct index_iterator<SparseMIABase<Derived> >: public index_iterator<Derived> {};

template<class Derived>
struct storage_iterator<SparseMIABase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<SparseMIABase<Derived> >: public const_storage_iterator<Derived> {};

template<class Derived>
struct Data<SparseMIABase<Derived> >: public Data<Derived> {};

template<class Derived>
struct Indices<SparseMIABase<Derived> >: public Indices<Derived> {};

}



//!  Base class for sparse multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  sparse multi-index array types. Provides operations and functions common to all sparse
  multi-index arrays.

  \tparam Derived   should only be a sparse multi-index array class.
*/
template <class Derived>
class SparseMIABase: public MIA<SparseMIABase<Derived > >
{
public:


    typedef typename internal::data_type<SparseMIABase>::type data_type;
    typedef typename internal::index_type<SparseMIABase>::type index_type;
    typedef typename internal::Data<SparseMIABase>::type Data;
    typedef typename internal::data_iterator<SparseMIABase>::type data_iterator;
    typedef typename internal::storage_iterator<SparseMIABase>::type storage_iterator;
    typedef typename internal::const_storage_iterator<SparseMIABase>::type const_storage_iterator;
    constexpr static size_t mOrder=internal::order<SparseMIABase>::value;
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
    SparseMIABase(Dims... dims): MIA<SparseMIABase<Derived > >(dims...) {}


    SparseMIABase(): MIA<SparseMIABase<Derived > >() {}


    SparseMIABase(const std::array<index_type,mOrder> &_dims): MIA<SparseMIABase<Derived > >(_dims) {}

//    template<class otherDerived>
//    bool operator==(const SparseMIABase<otherDerived>& otherMIA);
//
//    template<class otherDerived>
//    bool fuzzy_equals(const DenseMIABase<otherDerived> & otherMIA,data_type precision );
//
//    //! Returns scalar data at given indices
//    /*!
//        \param[in] indices variadic parameter. Will assert a compile error if size of indices!=mOrder or if Indices datatype are not convertible to index_type
//    */
//    template<typename... Indices>
//    const data_type& at(Indices... indices) const {
//        static_assert(internal::check_mia_constructor<SparseMIABase,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
//        std::array<index_type,internal::order<SparseMIABase>::value> temp = {{indices...}};
//        return (*(derived().data()))(temp);
//    }
//
//    //! Returns scalar data at given indices
//    /*!
//        \param[in] indices variadic parameter. Will assert a compile error if size of indices!=mOrder or if Indices datatype are not convertible to index_type
//    */
//    template<typename... Indices>
//    data_type& at(Indices... indices) {
//        static_assert(internal::check_mia_constructor<SparseMIABase,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
//        std::array<index_type,internal::order<SparseMIABase>::value> temp = {{indices...}};
//        return (*(derived().data()))(temp);
//    }
//
//    //!Assignment operator. Will call Derived's operator
//    template<class otherDerived>
//    SparseMIABase& operator=(const SparseMIABase<otherDerived>& otherMIA){
//        return derived()=otherMIA;
//    }
//
//
//    //!Assignment operator. Will call Derived's operator
//    SparseMIABase& operator=(const SparseMIABase& otherMIA){
//        return derived()=otherMIA.derived();
//    }
//
//    //!  Assignment based on given order. Will call Derived's operator
//    /*!
//
//        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
//        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up
//
//        \param[in] otherMIA the other MIA
//        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {3,1,2} this->at(1,2,3)==otherMIA.at(2,3,1).
//                                Will assert a compile failure is size of index_order is not the same as this->mOrder
//    */
//    template<class otherDerived,class index_param_type>
//    void assign(const MIA<otherDerived>& otherMIA,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order){
//        derived().assign(otherMIA.derived(),index_order);
//    }
//
//    //! Returns scalar data at given indices
//    const data_type& at(std::array<index_type, internal::order<SparseMIABase>::value> indices) const{
//
//        return (*(derived().data()))(indices);
//    }
//
//    //! Returns scalar data at given indices
//    data_type& at(std::array<index_type, internal::order<SparseMIABase>::value> indices){
//
//        return (*(derived().data()))(indices);
//    }
//
//    //! Returns scalar data at given linear index
//    const data_type& atIdx(index_type idx) const{
//
//        //return lin index
//        return *(derived().data_begin()+idx);
//    }
//
//    //! Returns scalar data at given linear index
//    data_type& atIdx(index_type idx) {
//
//        //return lin index
//        return *(derived().data_begin()+idx);
//    }

//    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom.
//    /*!
//        For DenseMIAs, this function calls toLatticeCopy
//    */
//    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
//    DenseLattice<data_type> toLatticeExpression(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const
//    {
//        return toLatticeCopy(row_indices, column_indices, tab_indices);
//
//    }
//
//    //! Flattens the MIA to a Lattice by creating a copy of the data.
//    /*!
//        \param[in] row_indices indices to map to the lattice rows - will perserve ordering
//        \param[in] column_indices indices to map to the lattice columns - will perserve ordering
//        \param[in] tab_indices indices to map to the lattice tabs - will perserve ordering
//        \return DenseLattice class that owns a copy of this's raw data
//    */
//    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
//    DenseLattice<data_type> toLatticeCopy(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices) const;
//
//    //! Performs destructive add (+=).
//    /*!
//        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is added like: this->at(x,y,z)+=b.at(y,z,x).
//        Will assert a compile failure is size of index_order is not the same as this->mOrder
//    */
//    template<class otherDerived,typename index_param_type>
//    DenseMIABase & plus_equal(const DenseMIABase<otherDerived> &b,const std::array<index_param_type,DenseMIABase::order>& index_order);
//
//
//    //! Performs destructive subtract (-=).
//    /*!
//        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is subtracted like: this->at(x,y,z)-=b.at(y,z,x).
//        Will assert a compile failure is size of index_order is not the same as this->mOrder
//    */
//    template<class otherDerived,typename index_param_type>
//    DenseMIABase & minus_equal(const DenseMIABase<otherDerived> &b,const std::array<index_param_type,DenseMIABase::order>& index_order);
//
//    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
//    template<typename otherDerived, typename Op,typename index_param_type>
//    void  merge(const DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,DenseMIABase::order>& index_order);

    data_iterator data_begin() const
    {


        return derived().data_begin();
    }

    data_iterator data_end() const
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



protected:




private:



};




/*! @} */

}

#endif // SPARSEMIABASE_H_INCLUDED
