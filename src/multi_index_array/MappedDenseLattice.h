// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


// Implementation of the DenseLattice class. Underlying data
// structure and algorithms are those provided by the
// EigenLib Library.

#ifndef MAPPEDDENSELATTICE_H
#define MAPPEDDENSELATTICE_H

#include <iostream>

#include <boost/shared_array.hpp>
#include <boost/multi_array.hpp>



#include "LibMIAUtil.h"
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

template<typename T>
struct const_data_iterator<MappedDenseLattice<T> >
{
    typedef const T* type;
};

template<typename T>
struct index_type<MappedDenseLattice<T> >
{
    typedef typename Eigen::Matrix<typename data_type<MappedDenseLattice<T>>::type,Eigen::Dynamic,Eigen::Dynamic>::Index type;
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
class MappedDenseLattice: public DenseLatticeBase<MappedDenseLattice<T> > 
{
public:
    //Lattice();
    //! Data type of individual elements
    typedef typename internal::data_type<MappedDenseLattice>::type data_type;
	typedef typename internal::index_type<MappedDenseLattice>::type index_type;
    typedef typename internal::Data<MappedDenseLattice>::type Data;
    typedef typename internal::storage_iterator<MappedDenseLattice>::type storage_iterator;
    typedef typename internal::const_storage_iterator<MappedDenseLattice>::type const_storage_iterator;
    typedef typename internal::data_iterator<MappedDenseLattice>::type data_iterator;
    typedef typename internal::const_data_iterator<MappedDenseLattice>::type const_data_iterator;
    typedef T* raw_pointer;



public:



    //!  Constructs MappedDenseLattice given already allocated and assigned memory.
    /*!

    \param[in] rawdata Pointer of already allocated memory. Internal pointer of class is set to rawdata.

    */
	MappedDenseLattice(T*rawdata, index_type _height, index_type _width, index_type _depth) : m_raw_data(rawdata), m_data(rawdata, boost::extents[_height][_width][_depth], boost::fortran_storage_order())
    {



        this->init( _height,  _width,  _depth);

    }

    template<class otherDerived>
    MappedDenseLattice& operator=(const DenseLatticeBase<otherDerived>& other);



	const data_type& operator()(index_type i, index_type j, index_type k) const
    {

        return m_data[i][j][k];
    }
	data_type& operator()(index_type i, index_type j, index_type k)
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

    const_data_iterator data_begin() const
    {
        return m_data.data();

    }

    const_data_iterator data_end() const
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

    //!  Releases the internal scalar data memory pointer.
    /*!
    Since this class doesn't own data, no other action is performed

    */
    raw_pointer release_memptr() const
    {


        return m_raw_data;


    }

	template<typename otherDerived1, typename otherDerived2>
	void perform_mult(const DenseLatticeBase<otherDerived1> & restrict_libmia a, const DenseLatticeBase<otherDerived2> & restrict_libmia b) restrict_libmia{

		if (a.depth() != b.depth())
			throw LatticeParameterException("Lattice depths must the same to be multipled.");
		if (a.width() != b.height())
			throw LatticeParameterException("Lattice width must be the same as Lattice height to be multiplied.");
		if (a.depth() != this->depth())
			throw LatticeParameterException("LHS Lattice depth must be the same as RHS Lattice depth.");
		if (b.width() != this->width())
			throw LatticeParameterException("LHS Lattice width must be the same as RHS Lattice width.");
		if (a.height() != this->height())
			throw LatticeParameterException("LHS Lattice height must be the same as RHS Lattice height.");

		for (int i = 0; i<this->depth(); i++)
		{

			this->tab_matrix(i).noalias() = (a.tab_matrix(i))*(b.tab_matrix(i));


		}
	}


protected:

private:

    //raw_pointer m_raw_data;
    raw_pointer m_raw_data;
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
