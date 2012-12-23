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

#ifndef MAPPEDDENSELATTICE_H
#define MAPPEDDENSELATTICE_H

#include <iostream>

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
struct data_type<MappedDenseLattice<T> >
{
    typedef T type;
};

template<typename T>
struct Data<MappedDenseLattice<T> >
{
    typedef typename boost::multi_array_ref<T,3> type;
};

template<typename T>
struct storage_iterator<MappedDenseLattice<T> >
{
    typedef typename Data<MappedDenseLattice<T> >::type::iterator type;
};

template<typename T>
struct const_storage_iterator<MappedDenseLattice<T> >
{
    typedef typename Data<MappedDenseLattice<T> >::type::const_iterator type;
};

template<typename T>
struct data_iterator<MappedDenseLattice<T> >
{
    typedef T* type;
};


}


/** \addtogroup lattice Lattice Classes
 *  @{
 */

//!  Lattice class for dense data. Maps already allocated data to a denselattice.
/*!
  Supports addition, multiplication, and solution of, possible over-determined, systems of
  linear equations. Class does not own internal scalar data - caller is responsible for
  freeing scalar data.

  \tparam T   the datatype of individual elements.
*/
template <class T>
class MappedDenseLattice: public DenseLatticeBase<MappedDenseLattice<T> > //, boost::multipliable<DenseLattice<T> >
{
public:
    //Lattice();
    //! Data type of individual elements
    typedef typename internal::data_type<MappedDenseLattice>::type data_type;
    typedef typename internal::Data<MappedDenseLattice>::type Data;
    typedef typename internal::storage_iterator<MappedDenseLattice>::type storage_iterator;
    typedef typename internal::const_storage_iterator<MappedDenseLattice>::type const_storage_iterator;
    typedef typename internal::data_iterator<MappedDenseLattice>::type data_iterator;
    typedef T* raw_pointer;



public:



    //!  Constructs MappedDenseLattice given already allocated and assigned memory.
    /*!

    \param[in] rawdata Pointer of already allocated memory. Internal pointer of class is set to rawdata.

    */
    MappedDenseLattice(T*rawdata,int _height, int _width, int _depth):  m_data(rawdata,boost::extents[_height][_width][_depth],boost::fortran_storage_order())
    {



        this->init( _height,  _width,  _depth);

    }

    template<class otherDerived>
    MappedDenseLattice& operator=(const DenseLatticeBase<otherDerived>& other);



    const data_type& operator()(int i, int j, int k) const
    {

        return m_data[i][j][k];
    }
    data_type& operator()(int i, int j, int k)
    {

        return m_data[i][j][k];
    }


    data_iterator data_begin()
    {
        return m_data.data();

    }

    data_iterator data_end()
    {
        return m_data.data()+size();

    }

    storage_iterator begin()
    {
        return m_data.begin();

    }

    storage_iterator end()
    {
        return m_data.end();

    }

    const_storage_iterator begin() const
    {
        return m_data.begin();

    }

    const_storage_iterator end() const
    {
        return m_data.end();

    }


    std::size_t size() const
    {

        return m_data.num_elements();

    }





protected:

private:

    //raw_pointer m_raw_data;

    Data m_data;




};




//****Operators*****

//self assignment is benign - no check is made
template <class T>
template<class otherDerived>
inline MappedDenseLattice<T>& MappedDenseLattice<T>::operator=(const DenseLatticeBase<otherDerived> &other)
{

    this->check_merge_dims(other);
    typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;

    for (int i=0;i<this->height();i++)
        for (int j=0;j<this->width();j++)
            for (int k=0;k<this->depth();k++)
                (*this)(i,j,k)= to_mdata_type::convert(other.derived()(i,j,k));



    return *this;



}


/*! @} */

}
typedef  LibMIA::MappedDenseLattice<double> dMapLattice;
typedef  LibMIA::MappedDenseLattice<int> iMapLattice;


#endif // MAPPEDDENSELATTICE_H
