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


#include "Util.h"
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





protected:

    typedef typename internal::data_type<DenseMIA>::type data_type;
    typedef typename internal::index_type<DenseMIA>::type index_type;
    typedef typename internal::Data<DenseMIA>::type Data;
    typedef typename internal::storage_iterator<DenseMIA>::type storage_iterator;
    typedef typename internal::const_storage_iterator<DenseMIA>::type const_storage_iterator;
    typedef typename internal::data_iterator<DenseMIA>::type data_iterator;
    typedef T* raw_pointer;
    constexpr static size_t value=_order;


private:
    typedef boost::shared_array<T> smart_pointer;
    smart_pointer m_smart_ptr;
    Data m_data;



public:

    struct null_deleter
    {
        void operator()(void const *) const
        {
        }
    };


    //!  Constructs DenseMIA of specified size.
    /*!

    */
    DenseMIA(std::array<index_type,_order> &_dims,T* scalar_data):DenseMIABase<DenseMIA<T,_order> >(_dims),m_smart_ptr(scalar_data),m_data(m_smart_ptr.get(),this->m_dims,boost::fortran_storage_order())
    {
    }

    //!  Constructs DenseMIA of specified size.
    /*!
        Extra template parameter is used to only enable constructor for

    */
    template<typename... Dims>
    DenseMIA(Dims... dims):DenseMIABase<DenseMIA<T,_order> >{dims...}, m_smart_ptr(new T[this->m_dimensionality]),m_data(m_smart_ptr.get(),this->m_dims,boost::fortran_storage_order())
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::cout<<this->m_dimensionality;



        //this->init( _height,  _width,  _depth);
    }

    //!  Constructs DenseMIA of specified size.
    /*!
        Initializes data to scalar array. Its size must be the same as m_dimensionality. Takes ownership of scalar_data

    */
    template<typename... Dims>
    DenseMIA(Dims... dims,T* scalar_data):DenseMIABase<DenseMIA<T,_order> >{dims...}, m_smart_ptr(scalar_data),m_data(m_smart_ptr.get(),this->m_dims,boost::fortran_storage_order())
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");



    }

    template<class...Ts>
    inline data_type at(Ts...ts){
        static_assert(internal::check_mia_constructor<DenseMIA,Ts...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::array<index_type,_order> temp ={{ts...}};
        return m_data(temp);


    }

    data_type atIdx(index_type idx) const{

        //return lin index
        return *(m_data.data()+idx);
    }

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
    data_iterator data_begin()
    {
        return m_data.data();

    }

    data_iterator data_end()
    {
        return m_data.data()+size();

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

        return m_data.num_elements();

    }



protected:

private:






};




////****Operators*****
//
////self assignment is benign - no check is made
//template <class T>
//template <class otherDerived>
//inline DenseMIA<T>& DenseMIA<T>::operator=(const DenseMIABase<otherDerived> &other)
//{
//
//    m_data.resize(boost::extents[other.height()][other.width()][other.depth()]);
//    this->init( other.height(), other.width(), other.depth());
//    typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;
//
//    for (int i=0;i<this->height();i++)
//        for (int j=0;j<this->width();j++)
//            for (int k=0;k<this->depth();k++)
//                (*this)(i,j,k)= to_mdata_type::convert(other.derived()(i,j,k));
//    return *this;
//
//
//
//}


/*! @} */

}

#endif // DENSEMIA_H_INCLUDED
