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




#include "LibMIAException.h"
#include "LibMIAUtil.h"
#include "FunctionUtil.h"
#include "libdivide.h"
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
	typedef typename std::make_unsigned<index_type>::type unsigned_index_type;
	
	typedef libdivide::divider<unsigned_index_type> fast_divisor;
	typedef std::array<unsigned_index_type, 3> accumulator_type;
	typedef std::array<unsigned_index_type, 3> multiplier_type;
	typedef std::array<fast_divisor, 3> fast_accumulator_type;


    Lattice() {}

    ~Lattice() {}


    Derived& derived() { return *static_cast<Derived*>(this); }
    /** \returns a const reference to the derived object */
    const Derived& derived() const { return *static_cast<const Derived*>(this); }


    inline index_type depth() const
    {
        return m_depth;
    }


    inline index_type width() const
    {
        return m_width;
    }
    inline index_type height() const
    {
        return m_height;
    }
    std::array<index_type,3> dims() const
    {
        std::array<index_type,3> ret={{m_height,m_width,m_depth}};
        return ret;
    }

    void set_height(index_type _height){
        m_height=_height;
		m_fast_height = fast_divisor(_height);
		m_tab_divisor = fast_divisor(this->width()*_height);
    }

    void set_width(index_type _width){
        m_width=_width;
		m_fast_width = fast_divisor(_width);
		m_tab_divisor = fast_divisor(this->height()*_width);
    }

    void set_depth(index_type _depth){
        m_depth=_depth;
		
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

    index_type sub2ind(index_type _row, index_type _column, index_type _tab)const{
        return _row+this->height()*(_column+this->width()*_tab);
    }


    /** Sets all lattice data to one.*/
    void ones(){
		std::fill(derived().data_begin(), derived().data_end(), data_type(1));
    }

    /** Sets all lattice data to one.*/
    void zeros(){
        std::fill ( derived().data_begin(), derived().data_end(), data_type(0));
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
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(LibMIA_gen(), uni_dist);
        typedef boost::numeric::converter<data_type,boost::uniform_real<>::result_type> to_mdata_type;
        for (data_iterator i=derived().data_begin();i<derived().data_end();i++){
            *i=to_mdata_type::convert(uni());
        }

    }
    SolveInfo mSolveInfo=NoInfo;

    SolveInfo solveInfo() const{
        return mSolveInfo;
    }
    void setSolveInfo(SolveInfo _solveInfo){
        mSolveInfo=_solveInfo;
    }

protected:

    /** Throws a LatticeParameterException if lattice multiplication operands have invalid dimensions.*/
    template<class otherDerived>
    void check_mult_dims(const Lattice<otherDerived> &rhs) const{
#ifdef LIBMIA_CHECK_DIMS
        if (rhs.derived().depth()!=derived().depth())
            throw LatticeParameterException("Lattice depths must the same to be multipled.");
        if (derived().width()!=rhs.derived().height())
            throw LatticeParameterException("LHS Lattice width must be the same as RHS Lattice height.");
#endif

    }

    /** Throws a LatticeParameterException if lattice solution operands have invalid dimensions.*/
    template<class otherDerived>
    void check_solve_dims(const Lattice<otherDerived> &rhs) const{
#ifdef LIBMIA_CHECK_DIMS
		if (rhs.derived().depth()!=derived().depth())
            throw LatticeParameterException("Lattice depths must be the same to solve systems of equations.");
        if (derived().height()!=rhs.derived().height())
            throw LatticeParameterException("LHS Lattice height must be the same as RHS Lattice height to solve system of equations.");
#endif

    }
    /** Throws a LatticeParameterException if lattice merging operands have invalid dimensions.*/
    template<class otherDerived>
    void check_merge_dims(const Lattice<otherDerived> &rhs) const{
#ifdef LIBMIA_CHECK_DIMS
		if (rhs.derived().depth()!=derived().depth() || rhs.derived().height()!=derived().height() || rhs.derived().width()!=derived().width())
            throw LatticeParameterException("Lattice dimensions must be identical for merger operation (+,-, etc).");
#endif

    }


    void init(index_type _height, index_type _width, index_type _depth, SolveInfo _solveInfo=NoInfo){
        m_height=_height;
        m_width=_width;
        m_depth=_depth;
		m_fast_height = fast_divisor(_height);
		m_fast_width = fast_divisor(_width);
		m_tab_divisor = fast_divisor(_height*_width);
        mSolveInfo=_solveInfo;

    }

	


	inline const fast_divisor& width_divisor() const
	{
		return m_fast_width;
	}
	inline const fast_divisor& height_divisor() const
	{
		return m_fast_height;
	}

	//! to get tab quickly
	inline const fast_divisor& tab_divisor() const
	{
		return m_tab_divisor;
	}


private:

    index_type m_height=0;
    index_type m_width=0;
    index_type m_depth=0;
	fast_divisor m_fast_height;
	fast_divisor m_fast_width;
	fast_divisor m_tab_divisor;
	




};


/*! @} */







} //namespace LibMIA
#endif // LATTICE_H
