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


}


/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */

//!  MIA class for dense data.
/*!
  Supports addition, multiplication, and solution of, possible over-determined, systems of
  linear equations. MIA owns underlying raw data, which can be reallocated and resized.

  \tparam T   the datatype of individual elements.
*/
template <class T, size_t _order>
class DenseMIA: public DenseMIABase<DenseMIA<T,_order> > //, boost::multipliable<DenseMIA<T> >
{





public:

    typedef typename internal::data_type<DenseMIA>::type data_type;
    typedef typename internal::index_type<DenseMIA>::type index_type;
    typedef typename internal::Data<DenseMIA>::type Data;
    typedef typename internal::storage_iterator<DenseMIA>::type storage_iterator;
    typedef typename internal::const_storage_iterator<DenseMIA>::type const_storage_iterator;
    typedef typename internal::data_iterator<DenseMIA>::type data_iterator;
    typedef T* raw_pointer;
    constexpr static size_t value=_order;


private:
    typedef std::unique_ptr<T []> smart_raw_pointer;
    typedef std::unique_ptr<Data> smart_data_pointer;
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

    DenseMIA():DenseMIABase<DenseMIA<T,_order> >(0),m_smart_raw_ptr(nullptr),m_Data(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order())),hasOwnership(true)
    {
    }

    //!  Constructs DenseMIA of specified size.
    /*!

    */
    template<class array_index_type>
    DenseMIA(const std::array<array_index_type,_order> &_dims,T* scalar_data,bool _ownership=true):DenseMIABase<DenseMIA<T,_order> >(_dims),hasOwnership(_ownership)
    {
        if(hasOwnership){
            m_smart_raw_ptr.reset(scalar_data);
            m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else{
            m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
            m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
            std::copy(scalar_data,scalar_data+this->m_dimensionality,this->data_begin());

        }

    }



    //!  Constructs DenseMIA of specified size.
    /*!
        Extra template parameter is used to only enable constructor for

    */
    template<typename... Dims>
    DenseMIA(Dims... dims):DenseMIABase<DenseMIA<T,_order> > {dims...}, m_smart_raw_ptr(new T[this->m_dimensionality]),m_Data(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order())),hasOwnership(true)
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        this->zeros();



        //this->init( _height,  _width,  _depth);
    }

    //!  Constructs DenseMIA of specified size.
    /*!
        Initializes data to scalar array. Its size must be the same as m_dimensionality. Takes ownership of scalar_data

    */
    template<typename... Dims>
    DenseMIA(T* scalar_data,bool _ownership,Dims... dims):DenseMIABase<DenseMIA<T,_order> > {dims...},hasOwnership(_ownership)
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        if(hasOwnership){
            m_smart_raw_ptr.reset(scalar_data);
            m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else{
            m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
            m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
            std::copy(scalar_data,scalar_data+this->m_dimensionality(),this->data_begin());

        }


    }

    template<typename otherDerived,typename index_param_type>
    void assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);

    template<typename otherDerived>
    DenseMIA& operator=(const DenseMIABase<otherDerived>& otherMIA);


    DenseMIA& operator=(const DenseMIA& otherMIA){
        return *this=static_cast<const DenseMIABase<DenseMIA>&>(otherMIA);
    }

    const smart_data_pointer & data() const{
        return m_Data;
    }



    T* raw_data_ptr() const{
        return m_smart_raw_ptr.get();
    }

    T* release_raw_data() {
        hasOwnership=false;
        return raw_data_ptr();
    }

//    template<class...Ts>
//    inline data_type at(Ts...ts)
//    {
//        static_assert(internal::check_mia_constructor<DenseMIA,Ts...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
//        std::array<index_type,_order> temp = {{ts...}};
//        return (*m_Data)(temp);
//
//
//    }
//
//    data_type atIdx(index_type idx) const
//    {
//
//        //return lin index
//        return *((*m_Data).data()+idx);
//    }

//    //!  Constructs empty DenseMIA.
//    //DenseMIA(): m_data() {}
//
//    //!  Constructs DenseMIA given already allocated and assigned memory.
//    /*!
//    Data is copied from raw pointer.
//    \param[in] rawdata Pointer of memory. Must have been allocated using new[].
//
//    */
//    DenseMIA(T*rawdata,int _height, int _width, int _depth):  m_smart_ptr(new T[_height*_width*_depth]),m_data(m_smart_ptr.get(),boost::extents[_height][_width][_depth],boost::fortran_storage_order())
//    {
//
//
//        std::copy(rawdata,rawdata+_height*_width*_depth,begin());
//        this->init( _height,  _width,  _depth);
//
//    }
//
//    //!  Copy constructor.
//    DenseMIA(const DenseMIA & other): m_smart_ptr(new T[other.height()*other.width()*other.depth()]), m_data(m_smart_ptr.get(),boost::extents[other.height()][other.width()][other.depth()],boost::fortran_storage_order())
//    {
//
//
//        m_data=other.m_data;
//        this->init( other.height(), other.width(), other.depth());
//
//    }
//
//    //!  Copy constructor.
//    template<class otherDerived>
//    DenseMIA(const DenseMIABase<otherDerived> & other): m_smart_ptr(new T[other.height()*other.width()*other.depth()]), m_data(m_smart_ptr.get(),boost::extents[other.height()][other.width()][other.depth()],boost::fortran_storage_order())
//    {
//
//        typedef  typename internal::data_type<otherDerived>::type other_data_type;
//        typedef boost::numeric::converter<data_type,other_data_type> to_mdata_type;
//        this->init( other.height(), other.width(), other.depth());
//        for (int i=0; i<this->height(); i++)
//            for(int j=0; j<this->width(); j++)
//                for(int k=0; k<this->depth(); k++)
//                    m_data[i][j][k]=to_mdata_type::convert(other.derived()(i,j,k));
//
//    }
//
//    //!  Releases the internal scalar data memory pointer.
//    /*!
//    After release, lattice is set to 0x0x0 and no longer refers to valid memory.
//    \return Pointer to scalar memory. Must be deallocated using delete[]. Caller is responsible for deallocation.
//
//    */
//    T* release_memptr(){
//        smart_pointer temp_ptr(0,null_deleter());
//        m_smart_ptr.swap(temp_ptr);
//        this->init(0,0,0);
//        return temp_ptr.get();
//
//
//    }
//
//    template<class otherDerived>
//    DenseMIA& operator=(const DenseMIABase<otherDerived>& b);
//
//
//
//    const data_type& operator()(int i, int j, int k) const
//    {
//
//        return m_data[i][j][k];
//    }
//    data_type& operator()(int i, int j, int k)
//    {
//
//        return m_data[i][j][k];
//    }
//
//
    data_iterator data_begin() const
    {
        return (*m_Data).data();

    }

    data_iterator data_end() const
    {
        return (*m_Data).data()+size();

    }
//
//    storage_iterator begin()
//    {
//        return m_data.begin();
//
//    }
//
//    storage_iterator end()
//    {
//        return m_data.end();
//
//    }
//
//    const_storage_iterator begin() const
//    {
//        return m_data.begin();
//
//    }
//
//    const_storage_iterator end() const
//    {
//        return m_data.end();
//
//    }
//
//
    //in dense case this is equivalent to dimensionality, but not in sparse case
    std::size_t size() const
    {

        return (*m_Data).num_elements();

    }

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
            m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
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

//template<class T, size_t _order>
//DenseMIA<T,_order>& DenseMIA<T,_order>::operator=(const DenseMIA<T,_order>& otherMIA)
//{
//    if(this->m_dims!=otherMIA.dims())
//        this->m_dims=otherMIA.dims();
//
//    if(this->m_dimensionality!=otherMIA.dimensionality()){
//
//        this->m_dimensionality=this->compute_dimensionality();
//
//        smart_raw_pointer temp_ptr(new T[this->m_dimensionality]);
//        m_smart_raw_ptr.swap(temp_ptr);
//        temp_ptr.release();
//        m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
//    }
//
//
//    for(auto it1=this->data_begin(),it2=otherMIA.data_begin();it1<this->data_end();++it1,++it2)
//        *it1=*it2;
//
//    return *this;
//
//}

template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void DenseMIA<T,_order>::assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");


    if(this->m_dimensionality!=otherMIA.dimensionality()){
        if(hasOwnership){
            internal::collect_dimensions_to_order(otherMIA,index_order,this->m_dims);
            this->m_dimensionality=otherMIA.dimensionality();
            smart_raw_pointer temp_ptr(new T[this->m_dimensionality]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();
            m_Data.reset(new Data(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    else
        internal::collect_dimensions_to_order(otherMIA,index_order,this->m_dims);

    typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;
    index_type curIdx=0;


    auto this_it=this->data_begin();
    for(auto it=otherMIA.data_begin(); it<otherMIA.data_end(); ++it)
    {
        *(this_it+sub2ind(ind2sub(curIdx++, otherMIA.dims()),index_order,this->dims()))=to_mdata_type::convert(*it);
    }


}




/*! @} */

}

#endif // DENSEMIA_H_INCLUDED
