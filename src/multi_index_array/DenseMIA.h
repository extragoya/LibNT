// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef DENSEMIA_H_INCLUDED
#define DENSEMIA_H_INCLUDED

#include <type_traits>
#include <iostream>
#include <algorithm>
#include <boost/shared_array.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>



#include "LibMiaException.h"
#include "Util.h"
#include "IndexUtil.h"
#include "DenseMIABase.h"



//\defgroup
namespace LibMIA
{

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */
namespace internal
{

template<typename T,size_t _order>
struct data_type<DenseMIA<T,_order> >
{
    typedef T type;
};

template<typename T,size_t _order>
struct index_type<DenseMIA<T,_order> >
{
    typedef long long type;
};

template<typename T,size_t _order>
struct order<DenseMIA<T,_order> >
{
    constexpr static size_t value=_order;
};


template<typename T,size_t _order>
struct Data<DenseMIA<T,_order> >
{
    typedef typename boost::multi_array_ref<T,_order> type;
};

template<typename T,size_t _order>
struct storage_iterator<DenseMIA<T,_order> >
{
    typedef typename Data<DenseMIA<T,_order> >::type::iterator type;
};

template<typename T,size_t _order>
struct const_storage_iterator<DenseMIA<T,_order> >
{
    typedef typename Data<DenseMIA<T,_order> >::type::const_iterator type;
};

template<typename T,size_t _order>
struct data_iterator<DenseMIA<T,_order> >
{
    typedef T* type;
};

template<typename T,size_t _order>
struct const_data_iterator<DenseMIA<T,_order> >
{
    typedef const T* type;
};

}




//!  MIA class for dense data.
/*!
  Supports addition, multiplication, and solution of, possibly over-determined, systems of
  linear equations. By default MIA will own underlying raw data, which can be reallocated and resized.

  \tparam T   the datatype of individual elements.
  \tparam _order   the order (number of indices) of the MIA.
*/
template <class T, size_t _order>
class DenseMIA: public DenseMIABase<DenseMIA<T,_order> >
{





public:

    //! raw data_type
    typedef typename internal::data_type<DenseMIA>::type data_type;
    //! raw index_type
    typedef typename internal::index_type<DenseMIA>::type index_type;
    //! data container type (that manages raw data)
    typedef typename internal::Data<DenseMIA>::type data_container_type;

    typedef typename internal::storage_iterator<DenseMIA>::type storage_iterator;
    typedef typename internal::const_storage_iterator<DenseMIA>::type const_storage_iterator;
    //! iterator type for iterating directly through raw data
    typedef typename internal::data_iterator<DenseMIA>::type data_iterator;
    //! raw data pointer type
    typedef T* raw_pointer;
    //! order of the MIA
    constexpr static size_t mOrder=_order;
    //! smart pointer type used to reference raw data
    typedef std::unique_ptr<T []> smart_raw_pointer;
    //! smart pointer type used to reference raw data container
    typedef std::unique_ptr<data_container_type> smart_data_pointer;

private:

    smart_raw_pointer m_smart_raw_ptr;
    smart_data_pointer m_Data;
    bool hasOwnership;


public:



    struct null_deleter
    {
        void operator()(void const *) const
        {
        }
    };

    //!  Constructs empty DenseMIA
    DenseMIA():DenseMIABase<DenseMIA<T,_order> >(),m_smart_raw_ptr(nullptr),m_Data(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order())),hasOwnership(true)
    {
    }

    //!  Constructs DenseMIA of specified size with a given raw data pointer.
    /*!
        \tparam[in] _dims The dimensions size of data. Will assert a compile failure is size is different than _order
        \tparam[in] scalar_data Raw pointer of scalar data
        \tparam[in] _ownership If false, DenseMIA will not own data and caller is responsible for deleting scalar_data. Also DenseMIA will not be
                                allowed to resize scalar_data, but can assign individual entries.
    */
    template<class array_index_type>
    DenseMIA(const std::array<array_index_type,_order> &_dims,T* scalar_data,bool _ownership=true):DenseMIABase<DenseMIA<T,_order> >(_dims),hasOwnership(_ownership)
    {
        if(hasOwnership){
            m_smart_raw_ptr.reset(scalar_data);
            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else{
            m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
            std::copy(scalar_data,scalar_data+this->m_dimensionality,this->data_begin());

        }

    }

    //!  Copy constructor.
    /*!


    */
    DenseMIA(const DenseMIA& otherMIA):DenseMIABase<DenseMIA<T,_order> >(otherMIA.dims())
    {


        m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
        m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        std::copy(otherMIA.data_begin(),otherMIA.data_end(),this->data_begin());


    }

    //!  Copy constructor.
    /*!
        If otherMIA's datatype is different than this->data_type, then individual entries will be converted.

    */
    template<class otherDerived>
    DenseMIA(const DenseMIABase<otherDerived>& otherMIA):DenseMIABase<DenseMIA<T,_order> >(otherMIA.dims())
    {



        static_assert(internal::order<otherDerived>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");
        m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
        m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        size_t idx=0;
        std::for_each( otherMIA.data_begin(), otherMIA.data_end(), [this,&idx] (typename otherDerived::data_type val)
        {
            *(this->data_begin()+idx++)=this->convert<typename otherDerived::data_type>(val);
        } );




    }


    //!  Constructs DenseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<typename... Dims>
    DenseMIA(Dims... dims):DenseMIABase<DenseMIA<T,_order> > {dims...}, m_smart_raw_ptr(new T[this->m_dimensionality]),m_Data(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order())),hasOwnership(true)
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        this->zeros();

    }

    //!  Constructs DenseMIA of specified size with a given raw data pointer.
    /*!

        \param[in] scalar_data Raw pointer of scalar data
        \param[in] _ownership If false, DenseMIA will not own data and caller is responsible for deleting scalar_data. Also DenseMIA will not be
                                allowed to resize scalar_data, but can assign individual entries.
        \param[in] _dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than mOrder or if datatype of variables making
                            up _dims are not convertible to index_type
    */
    template<typename... Dims>
    DenseMIA(T* scalar_data,bool _ownership,Dims... dims):DenseMIABase<DenseMIA<T,_order> > {dims...},hasOwnership(_ownership)
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        if(hasOwnership){
            m_smart_raw_ptr.reset(scalar_data);
            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else{
            m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
            std::copy(scalar_data,scalar_data+this->m_dimensionality(),this->data_begin());

        }


    }

    //!  Assignment based on given order.
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up

        \param[in] otherMIA the other MIA
        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {2,0,1} this->at(x,y,z)==otherMIA.at(y,z,x).

    */
    template<typename otherDerived,typename index_param_type>
    void assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);

    //!  Assignment operator.
    /*!
        \param[in] otherMIA
        If the data_type of otherMIA is not the same as this, the scalar data will be converted.
    */
    template<typename otherDerived>
    DenseMIA& operator=(const DenseMIABase<otherDerived>& otherMIA);

    //!  Straight-out assignment.
    DenseMIA& operator=(const DenseMIA& otherMIA){
        return *this=static_cast<const DenseMIABase<DenseMIA>&>(otherMIA);
    }

    //! Returns a smart pointer to the data container used
    const smart_data_pointer & data() const{
        return m_Data;
    }


    //! Returns a raw pointer to the scalar data
    T* raw_data_ptr() const{
        return m_smart_raw_ptr.get();
    }

    //! Returns a raw pointer to the scalar data and releases ownership - caller must deallocate data using delete[]
    T* release_raw_data() {
        hasOwnership=false;
        return raw_data_ptr();
    }

    //! Iterator to the beginning of the raw data
    data_iterator data_begin() const
    {
        return (*m_Data).data();

    }

    //! Iterator to the end of the raw data
    data_iterator data_end() const
    {
        return (*m_Data).data()+size();

    }

    //! Converts a scalar value to data_type
    /*!
        \tparam from_data_type the data_type you are converting from
    */
    template<class from_data_type>
    data_type convert(const from_data_type from) const{
        using namespace boost::numeric;
        typedef converter<data_type,from_data_type,conversion_traits<data_type,from_data_type>,def_overflow_handler,RoundEven<from_data_type>> to_mdata_type;
        return to_mdata_type::convert(from);
    }

    //! Returns size of raw data. For dense cases, this is the same as dimensionality
    std::size_t size() const
    {

        return (*m_Data).num_elements();

    }



    //! If this has ownership, the raw data will be deallocated
    ~DenseMIA(){
        if(!hasOwnership)
            m_smart_raw_ptr.release();
    }

protected:

private:






};


template<class T, size_t _order>
template<typename otherDerived>
DenseMIA<T,_order>& DenseMIA<T,_order>::operator=(const DenseMIABase<otherDerived>& otherMIA)
{
    static_assert(_order==internal::order<otherDerived>::value,"Orders of MIAs must be the same to be assigned");


    if(this->m_dimensionality!=otherMIA.dimensionality()){
        if(hasOwnership){
            this->m_dimensionality=otherMIA.dimensionality();
            this->m_dims=otherMIA.dims();
            smart_raw_pointer temp_ptr(new T[this->m_dimensionality]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();
            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    else if(this->m_dims!=otherMIA.dims())
        this->m_dims=otherMIA.dims();

    typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;
    for(auto it1=this->data_begin(),it2=otherMIA.data_begin();it1<this->data_end();++it1,++it2)
        *it1=to_mdata_type::convert(*it2);

    return *this;

}






template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void DenseMIA<T,_order>::assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");


    if(this->m_dimensionality!=otherMIA.dimensionality()){
        if(hasOwnership){
            internal::reorder_from(otherMIA.dims(),index_order,this->m_dims);
            this->m_dimensionality=otherMIA.dimensionality();
            smart_raw_pointer temp_ptr(new T[this->m_dimensionality]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();
            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    else
        internal::reorder_from(otherMIA.dims(),index_order,this->m_dims);


    index_type curIdx=0;

    auto other_it=otherMIA.data_begin();
    for(auto this_it=this->data_begin(); this_it<this->data_end(); ++this_it)
    {
        *this_it=this->convert(*(other_it+internal::sub2ind(this->ind2sub(curIdx++),index_order,otherMIA.dims())));

    }


}




/*! @} */

}

#endif // DENSEMIA_H_INCLUDED
