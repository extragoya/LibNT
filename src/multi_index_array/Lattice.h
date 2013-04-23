// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.




#ifndef LATTICE_H
#define LATTICE_H

#include <stdlib.h>

#include <algorithm>

#include <boost/operators.hpp>


#include "LibMIAException.h"
#include "Util.h"

//, boost::multipliable<Derived>
//
namespace LibMIA
{

namespace internal{

template<class Derived>
struct index_type<Lattice<Derived> >: public index_type<Derived> {};

template<class Derived>
struct data_type<Lattice<Derived> >: public data_type<Derived> {};

template<class Derived>
struct data_iterator<Lattice<Derived> >: public data_iterator<Derived> {};

}

/** \addtogroup lattice Lattice Classes
*  @{
*/


//!  Base class for lattice classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  lattice types. Provides common operations and functions.

  \tparam Derived   should only be DenseLatticeBase or SparseLatticeBase type.
*/
template
<

    class Derived
>
class Lattice  //boost::additive<Derived>
{
public:




    typedef typename internal::data_type<Lattice>::type data_type;
    typedef typename internal::index_type<Lattice>::type index_type;
    typedef typename internal::data_iterator<Lattice>::type data_iterator;


    Lattice() {}

    ~Lattice() {}


    Derived& derived() { return *static_cast<Derived*>(this); }
    /** \returns a const reference to the derived object */
    const Derived& derived() const { return *static_cast<const Derived*>(this); }


    index_type depth() const
    {
        return m_depth;
    }
    index_type width() const
    {
        return m_width;
    }
    index_type height() const
    {
        return m_height;
    }

    index_type dimensionality() const
    {
        return m_height*m_width*m_depth;
    }

    //! Converts a scalar value to data_type
    /*!
        \tparam from_data_type the data_type you are converting from
    */
    template<class from_data_type>
    static data_type convert(const from_data_type from){

        return internal::convert<data_type,from_data_type>(from);
    }


    /** Sets all lattice data to one.*/
    void ones(){
        std::fill ( derived().data_begin(), derived().data_end(), 1);
    }

    /** Sets all lattice data to one.*/
    void zeros(){
        std::fill ( derived().data_begin(), derived().data_end(), 0);
    }

    //!  Sets lattice data to uniformly distributed random values.
    /*!
    If lattice is sparse, only elements already designated as nonzero are set to random values.
    Range is specified using parameters. Will throw a LatticeParameterException exception if \f$low>high\f$.

    */
    void randu(int low=0, int high=1){
        if (low>=high){
            throw LatticeParameterException("Lower bound of random numbers must be stricly smaller than upper bound.");
        }
        boost::uniform_real<> uni_dist(low,high);
        boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(gen, uni_dist);
        typedef boost::numeric::converter<data_type,boost::uniform_real<>::result_type> to_mdata_type;
        for (data_iterator i=derived().data_begin();i<derived().data_end();i++){
            *i=to_mdata_type::convert(uni());
        }

    }


protected:

    /** Throws a LatticeParameterException if lattice multiplication operands have invalid dimensions.*/
    template<class otherDerived>
    void check_mult_dims(const Lattice<otherDerived> &rhs) const{

        if (rhs.derived().depth()!=derived().depth())
            throw LatticeParameterException("Lattice depths must the same to be multipled.");
        if (derived().width()!=rhs.derived().height())
            throw LatticeParameterException("LHS Lattice height must be the same as RHS Lattice width.");

    }

    /** Throws a LatticeParameterException if lattice solution operands have invalid dimensions.*/
    template<class otherDerived>
    void check_solve_dims(const Lattice<otherDerived> &rhs) const{
        if (rhs.derived().depth()!=derived().depth())
            throw LatticeParameterException("Lattice depths must be the same to solve systems of equations.");
        if (derived().height()!=rhs.derived().height())
            throw LatticeParameterException("LHS Lattice height must be the same as RHS Lattice height to solve system of equations.");

    }
    /** Throws a LatticeParameterException if lattice merging operands have invalid dimensions.*/
    template<class otherDerived>
    void check_merge_dims(const Lattice<otherDerived> &rhs) const{
        if (rhs.derived().depth()!=derived().depth() || rhs.derived().height()!=derived().height() || rhs.derived().width()!=derived().width())
            throw LatticeParameterException("Lattice dimensions must be identical for merger operation (+,-, etc).");


    }


    void init(index_type _height, index_type _width, index_type _depth){
        m_height=_height;
        m_width=_width;
        m_depth=_depth;


    }

protected:


    index_type m_height;
    index_type m_width;
    index_type m_depth;





};


/*! @} */







} //namespace LibMIA
#endif // LATTICE_H
