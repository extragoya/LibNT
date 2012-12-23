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



#ifndef LATTICE_H
#define LATTICE_H

#include <stdlib.h>
#include <time.h>
#include <algorithm>

#include <boost/operators.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "LatticeException.h"
#include "Util.h"

//, boost::multipliable<Derived>
//
namespace LibMIA
{
boost::random::mt19937 gen(time(0));
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




    typedef typename internal::data_type<Derived>::type data_type;
    typedef typename internal::index_type<Derived>::type index_type;
    typedef typename internal::data_iterator<Derived>::type data_iterator;


    Lattice() {}

    ~Lattice() {}


    Derived& derived() { return *static_cast<Derived*>(this); }
    /** \returns a const reference to the derived object */
    const Derived& derived() const { return *static_cast<const Derived*>(this); }


    int depth() const
    {
        return m_depth;
    }
    int width() const
    {
        return m_width;
    }
    int height() const
    {
        return m_height;
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
        std::cout << rhs.derived().depth() << " " << derived().depth() << "\n";
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

private:


    index_type m_height;
    index_type m_width;
    index_type m_depth;





};


/*! @} */







} //namespace LibMIA
#endif // LATTICE_H
