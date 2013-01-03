// Copyright (C) 2011 Adam P. Harrison
// http://www.ualberta.ca/~apharris/
//
// This file is part of the MIA library, developed at
// the Electronic Imaging Lab, University of Alberta,
// Edmonton, Canada
// http://www.ece.ualberta.ca/~djoseph/
//
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)

// Implementation of the DenseLattice class. Underlying data
// structure and algorithms are those provided by the
// EigenLib Library.

#ifndef DENSELATTICE_H
#define DENSELATTICE_H

#include <iostream>
#include <memory>

#include <boost/shared_array.hpp>
#include <boost/multi_array.hpp>



#include "Util.h"
#include "DenseLatticeBase.h"


//\defgroup
namespace LibMIA
{


namespace internal
{

template<typename T>
struct data_type<DenseLattice<T> >
{
    typedef T type;
};



template<typename T>
struct Data<DenseLattice<T> >
{
    typedef typename boost::multi_array_ref<T,3> type;
};

template<typename T>
struct storage_iterator<DenseLattice<T> >
{
    typedef typename Data<DenseLattice<T> >::type::iterator type;
};

template<typename T>
struct const_storage_iterator<DenseLattice<T> >
{
    typedef typename Data<DenseLattice<T> >::type::const_iterator type;
};

template<typename T>
struct data_iterator<DenseLattice<T> >
{
    typedef T* type;
};


}


/** \addtogroup lattice Lattice Classes
 *  @{
 */

//!  Lattice class for dense data.
/*!
  Supports addition, multiplication, and solution of, possible over-determined, systems of
  linear equations. Lattice owns underlying raw data, which can be reallocated and resized.

  \tparam T   the datatype of individual elements.
*/
template <class T>
class DenseLattice: public DenseLatticeBase<DenseLattice<T> > //, boost::multipliable<DenseLattice<T> >
{
public:

    typedef typename internal::data_type<DenseLattice>::type data_type;
    typedef typename internal::Data<DenseLattice>::type Data;
    typedef typename internal::storage_iterator<DenseLattice>::type storage_iterator;
    typedef typename internal::const_storage_iterator<DenseLattice>::type const_storage_iterator;
    typedef typename internal::data_iterator<DenseLattice>::type data_iterator;
    typedef T* raw_pointer;

private:
    typedef std::unique_ptr<T []> smart_raw_pointer;
    typedef std::unique_ptr<Data> smart_data_pointer;

public:

    struct null_deleter
    {
        void operator()(void const *) const
        {
        }
    };

    //!  Constructs DenseLattice of specified size, initialized to all zeros.
    DenseLattice(int _height, int _width, int _depth) : m_smart_raw_ptr(new T[_height*_width*_depth]), m_Data(new Data(m_smart_raw_ptr.get(),boost::extents[_height][_width][_depth],boost::fortran_storage_order()))
    {

        this->zeros();
        this->init( _height,  _width,  _depth);
    }

    //!  Constructs empty DenseLattice.
    DenseLattice():  m_smart_raw_ptr(),m_Data(new Data(m_smart_raw_ptr.get(),boost::extents[0][0][0],boost::fortran_storage_order())) {}

    //!  Constructs DenseLattice given already allocated and assigned memory.
    /*!
    Data is copied from raw pointer.
    \param[in] rawdata Pointer of memory.

    */
    DenseLattice(T*rawdata,int _height, int _width, int _depth):  m_smart_raw_ptr(new T[_height*_width*_depth]),m_Data(new Data(m_smart_raw_ptr.get(),boost::extents[_height][_width][_depth],boost::fortran_storage_order()))
    {


        std::copy(rawdata,rawdata+_height*_width*_depth,begin());
        this->init( _height,  _width,  _depth);

    }

    //!  Copy constructor.
    DenseLattice(const DenseLattice & other): m_smart_raw_ptr(new T[other.height()*other.width()*other.depth()]), m_Data(new Data(m_smart_raw_ptr.get(),boost::extents[other.height()][other.width()][other.depth()],boost::fortran_storage_order()))
    {


        *m_Data=*(other.m_Data);
        this->init( other.height(), other.width(), other.depth());

    }

    //!  Move constructor.
    DenseLattice(DenseLattice && other) :m_smart_raw_ptr(),m_Data(new Data(m_smart_raw_ptr.get(),boost::extents[0][0][0],boost::fortran_storage_order()))
    {

        m_smart_raw_ptr.swap(other.m_smart_raw_ptr);
        m_Data.swap(other.m_Data);
        this->init( other.height(), other.width(), other.depth());

    }



    //TODO move constructor and move assignment

    //!  Copy constructor.
    template<class otherDerived>
    DenseLattice(const DenseLatticeBase<otherDerived> & other): m_smart_raw_ptr(new T[other.height()*other.width()*other.depth()]), m_Data(new Data(m_smart_raw_ptr.get(),boost::extents[other.height()][other.width()][other.depth()],boost::fortran_storage_order()))
    {

        typedef  typename internal::data_type<otherDerived>::type other_data_type;
        typedef boost::numeric::converter<data_type,other_data_type> to_mdata_type;
        this->init( other.height(), other.width(), other.depth());
        for (int i=0; i<this->height(); i++)
            for(int j=0; j<this->width(); j++)
                for(int k=0; k<this->depth(); k++)
                    *m_Data[i][j][k]=to_mdata_type::convert(other.derived()(i,j,k));

    }

    //!  Releases the internal scalar data memory pointer.
    /*!
    After release, lattice is set to 0x0x0 and no longer refers to valid memory.
    \return Pointer to scalar memory. Must be deallocated using delete[]. Caller is responsible for deallocation.

    */
    T* release_memptr()
    {

        smart_raw_pointer temp_ptr;
        m_smart_raw_ptr.swap(temp_ptr);
        m_Data.reset(new Data(m_smart_raw_ptr.get(),boost::extents[0][0][0],boost::fortran_storage_order()));
        return temp_ptr.release();


    }



    template<class otherDerived>
    DenseLattice& operator=(const DenseLatticeBase<otherDerived>& b);

    //Move assignment
    DenseLattice& operator=(DenseLattice && other)
    {
        if (this != &other)
        {
            m_smart_raw_ptr.swap(other.m_smart_raw_ptr);
            m_Data.swap(other.m_Data);
            this->init( other.height(), other.width(), other.depth());
        }
        return *this;

    }

    const data_type& operator()(int i, int j, int k) const
    {

        return (*m_Data)[i][j][k];
    }
    data_type& operator()(int i, int j, int k)
    {

        return (*m_Data)[i][j][k];
    }


    data_iterator data_begin()
    {
        return (*m_Data).data();

    }

    data_iterator data_end()
    {
        return (*m_Data).data()+size();

    }

    storage_iterator begin()
    {
        return (*m_Data).begin();

    }

    storage_iterator end()
    {
        return (*m_Data).end();

    }

    const_storage_iterator begin() const
    {
        return (*m_Data).begin();

    }

    const_storage_iterator end() const
    {
        return (*m_Data).end();

    }


    std::size_t size() const
    {

        return (*m_Data).num_elements();

    }



protected:

private:

    //raw_pointer m_raw_data;
    smart_raw_pointer m_smart_raw_ptr;
    smart_data_pointer m_Data;
    //Data m_data;




};




//****Operators*****


template <class T>
template <class otherDerived>
inline DenseLattice<T>& DenseLattice<T>::operator=(const DenseLatticeBase<otherDerived> &other)
{

    if (this != &other)
    {
        smart_raw_pointer temp_ptr(new T[other.height()*other.width()*other.depth()]);
        m_smart_raw_ptr.swap(temp_ptr);
        m_Data.reset(new Data(m_smart_raw_ptr.get(),boost::extents[other.height()][other.width()][other.depth()],boost::fortran_storage_order()));
        this->init( other.height(), other.width(), other.depth());
        typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;

        for (int i=0; i<this->height(); i++)
            for (int j=0; j<this->width(); j++)
                for (int k=0; k<this->depth(); k++)
                    (*this)(i,j,k)= to_mdata_type::convert(other.derived()(i,j,k));
    }
    return *this;



}


/*! @} */

}
typedef  LibMIA::DenseLattice<double> dLattice;
typedef  LibMIA::DenseLattice<int> iLattice;


#endif // DENSELATTICE_H
