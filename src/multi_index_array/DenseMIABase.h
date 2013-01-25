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


#include "Util.h"
#include "MIA.h"
#include "DenseLattice.h"



namespace LibMIA
{

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
struct storage_iterator<DenseMIABase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<DenseMIABase<Derived> >: public const_storage_iterator<Derived> {};

}

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */

//!  Base class for dense multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  dense multi-index array types. Provides operations and functions common to all dense
  multi-index array.

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
    typedef typename internal::storage_iterator<Derived>::type storage_iterator;
    typedef typename internal::const_storage_iterator<Derived>::type const_storage_iterator;

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
    DenseMIABase(Dims... dims): MIA<DenseMIABase<Derived > >(dims...) {}


    DenseMIABase(): MIA<DenseMIABase<Derived > >() {}


    DenseMIABase(std::array<index_type,internal::order<DenseMIABase>::value> &_dims): MIA<DenseMIABase<Derived > >(_dims) {}

    template<class otherDerived>
    bool operator==(const DenseMIABase<otherDerived>& otherMIA);

    template<typename... Indices>
    const data_type& at(Indices... indices) const {
        static_assert(internal::check_mia_constructor<DenseMIABase,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::array<index_type,internal::order<DenseMIABase>::value> temp = {{indices...}};
        return (*(derived().data()))(temp);
    }

    template<typename... Indices>
    data_type& at(Indices... indices) {
        static_assert(internal::check_mia_constructor<DenseMIABase,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::array<index_type,internal::order<DenseMIABase>::value> temp = {{indices...}};
        return (*(derived().data()))(temp);
    }

    template<class otherDerived>
    DenseMIABase& operator=(const DenseMIABase<otherDerived>& otherMIA){
        return derived()=otherMIA;
    }


    DenseMIABase& operator=(const DenseMIABase& otherMIA){
        return derived()=otherMIA.derived();
    }

    template<class otherDerived,class index_param_type>
    void assign(const MIA<otherDerived>& otherMIA,const std::array<index_param_type,internal::order<DenseMIABase>::value>& index_order){
        derived().assign(otherMIA.derived(),index_order);
    }

    const data_type& at(std::array<index_type, internal::order<DenseMIABase>::value> indices) const{

        return (*(derived().data()))(indices);
    }

    data_type& at(std::array<index_type, internal::order<DenseMIABase>::value> indices){

        return (*(derived().data()))(indices);
    }

    const data_type& atIdx(index_type idx) const{

        //return lin index
        return *(derived().data_begin()+idx);
    }

    data_type& atIdx(index_type idx) {

        //return lin index
        return *(derived().data_begin()+idx);
    }

    template<typename R,typename C, typename T>
    DenseLattice<data_type> toLatticeExpression(internal::sequence_array<R> row_indices, internal::sequence_array<C> column_indices,internal::sequence_array<T> tab_indices) const
    {
        return toLatticeCopy(row_indices, column_indices, tab_indices);

    }

    template<typename R,typename C, typename T>
    DenseLattice<data_type> toLatticeCopy(internal::sequence_array<R> row_indices, internal::sequence_array<C> column_indices,internal::sequence_array<T> tab_indices) const;


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




template<typename Derived>
template< typename R, typename C, typename T>
auto DenseMIABase<Derived>::toLatticeCopy(internal::sequence_array<R> row_indices, internal::sequence_array<C> column_indices,internal::sequence_array<T> tab_indices) const ->DenseLattice<data_type>
{

    //statically check number of indices match up
    size_t row_size=1, column_size=1, tab_size=1;
    std::array<index_type, boost::mpl::size<R>::type::value> row_dims;
    std::array<index_type, boost::mpl::size<C>::type::value> column_dims;
    std::array<index_type, boost::mpl::size<T>::type::value> tab_dims;
    size_t idx=0;
    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    for(auto _row: row_indices)
    {

        row_size*=this->m_dims[_row];
        row_dims[idx++]=this->m_dims[_row];
    }

    idx=0;

    for(auto _column: column_indices)
    {

        column_size*=this->m_dims[_column];
        column_dims[idx++]=this->m_dims[_column];
    }
    idx=0;

    for(auto _tab: tab_indices)
    {

        tab_size*=this->m_dims[_tab];
        tab_dims[idx++]=this->m_dims[_tab];
    }

    //std::cout<< "Tab dims " << tab_dims[0] << "\n";
    DenseLattice<data_type> lat(row_size, column_size, tab_size);



    index_type row_idx=0, column_idx=0,tab_idx=0;
    for (index_type k=0;k<tab_size;++k){
        tab_idx=sub2ind(ind2sub(k,tab_dims),tab_indices,this->m_dims);
        for (index_type j=0;j<column_size;++j){
            column_idx=tab_idx+sub2ind(ind2sub(j,column_dims),column_indices,this->m_dims);
            for (index_type i=0;i<row_size;++i){
                row_idx=column_idx+sub2ind(ind2sub(i,row_dims),row_indices,this->m_dims);

                lat(i,j,k)=this->atIdx(row_idx);

            }
        }


    }
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



/*! @} */

}






#endif // DENSEMIABASE_H_INCLUDED
