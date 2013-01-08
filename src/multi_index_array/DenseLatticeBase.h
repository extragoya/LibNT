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
// Armadillo C++ Library. Uses the boost::operators class to
// implement operators.

#ifndef DENSELATTICEBASE_H
#define DENSELATTICEBASE_H
#include <iostream>
#include <fstream>
#include <string>


#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/type_traits/is_floating_point.hpp>

#include <Eigen/Dense>

#include "Util.h"
#include "Lattice.h"
#include "Packer.h"



namespace LibMIA
{

namespace internal
{



template<class Derived>
struct data_type<DenseLatticeBase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<DenseLatticeBase<Derived> >
{
    typedef typename Eigen::Matrix<typename data_type<Derived>::type,Eigen::Dynamic,Eigen::Dynamic>::Index type;
};


template<class Derived>
struct data_iterator<DenseLatticeBase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct storage_iterator<DenseLatticeBase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<DenseLatticeBase<Derived> >: public const_storage_iterator<Derived> {};

}

/** \addtogroup lattice Lattice Classes
 *  @{
 */

//!  Base class for dense lattice classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  dense lattice types. Provides operations and functions common to all dense
  lattices.

  \tparam Derived   should only be a dense lattice class.
*/
template <class Derived>
class DenseLatticeBase: public Lattice<DenseLatticeBase<Derived > >//, boost::multipliable<DenseLattice<T> >
{
public:


    typedef typename internal::data_type<Derived>::type data_type;
    typedef typename internal::Data<Derived>::type Data;
    typedef typename internal::data_iterator<Derived>::type data_iterator;
    typedef typename internal::storage_iterator<Derived>::type storage_iterator;
    typedef typename internal::const_storage_iterator<Derived>::type const_storage_iterator;
    typedef typename Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> > matrix_type;
    typedef const typename Eigen::Map<const Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> > const_matrix_type;
    typedef typename Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, 1> > vector_type;
    typedef const typename Eigen::Map<const Eigen::Matrix<data_type, Eigen::Dynamic, 1> > const_vector_type;

    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }
    /** \returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }


    //!  Performs lattice product.
    /*!
        c=a*b. Performs matrix/matrix product across each of the operand types. Uses <a href="http://eigen.tuxfamily.org">Eigen</a> to perform
        matrix products. Operation is only defined for two dense lattice of the same scalar type (which may change as EigenLib evolves).

        \returns DenseLattice<data_type> regardless of the operand types.
    */
    template<class otherDerived>
    typename DenseProductReturnType<Derived,otherDerived>::type
    operator*(const DenseLatticeBase<otherDerived> &b) const;

    //!  Solves linear system of equations.
    /*!
    x=a.solve(b), solves a sequence of systems of linear equations represented by a*x=b. Uses <a href="http://eigen.tuxfamily.org">Eigen</a> to compute
        the solution. System can be over-determined. Operation is only defined for two dense lattice of the same scalar type (which may change as EigenLib evolves). Scalar type
        must be floating point.

        \returns DenseLattice<data_type> regardless of the operand types.
    */
    template<class otherDerived>
    typename DenseSolveReturnType<Derived,otherDerived>::type
    solve(const DenseLatticeBase<otherDerived> &b) const;

    //!  Similar to solve, but assumes non-singular 'a'.
    template<class otherDerived>
    void inverse_solve(const DenseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const;

    //!  Similar to solve, but assumes overdetermined 'a'.
    template<class otherDerived>
    void lsqr_solve(const DenseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const;

    void print() const;

    bool save(const std::string & _filename) const;

    template<class otherDerived>
    bool operator==(const DenseLatticeBase<otherDerived> &b) const;

    //!  Sets each tab to an identity matrix.
    /*!
        If non-square, it follows <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1MatrixBase.html#a0650b65c6ae6c3d19a138b72a6d68568">this</a> format.
    */
    void eye();

    //!  Destructive add.
    /*!
        \returns DenseLattice<data_type>. If data_type is a subtype of b's data_type, precision may be lost. Conversion executed using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    DenseLatticeBase & operator+=(const DenseLatticeBase<otherDerived> &b)
    {
        this->check_merge_dims(b);
        std::plus<data_type> op;
        merge(b,op);
        return *this;

    }

    //!  Destructive subtract.
    /*!
        \returns DenseLattice<data_type>. If data_type is a subtype of b's data_type, precision may be lost. Conversion executed using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    DenseLatticeBase & operator-=(const DenseLatticeBase<otherDerived> &b)
    {
        this->check_merge_dims(b);
        std::minus<data_type> op;
        merge(b,op);
        return *this;

    }

    //!  Binary addition.
    /*!
        \returns DenseLattice<super_data_type>. Super_data_type is determined using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    typename DenseMergeReturnType<Derived,otherDerived>::type operator+(DenseLatticeBase<otherDerived> &b)
    {
        typedef typename DenseMergeReturnType<Derived,otherDerived>::type RType;
        this->check_merge_dims(b);
        RType c(*this);
        c+=b;
        return c;

    }

    //!  Binary subtraction.
    /*!
        \returns DenseLattice<super_data_type>. Super_data_type is determined using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    typename DenseMergeReturnType<Derived,otherDerived>::type operator-(DenseLatticeBase<otherDerived> &b)
    {
        typedef typename DenseMergeReturnType<Derived,otherDerived>::type RType;
        this->check_merge_dims(b);
        RType c(*this);
        c-=b;
        return c;

    }



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

    //!  Wraps specified tab with an <a href="http://eigen.tuxfamily.org">Eigen</a> matrix.
    matrix_type tab_matrix(int _tab){
        return matrix_type(&(derived()(0,0,_tab)),this->height(),this->width());

    }

    const_matrix_type tab_matrix(int _tab) const{
        return const_matrix_type(&(derived()(0,0,_tab)),this->height(),this->width());

    }

    //!  Wraps specified column with an <a href="http://eigen.tuxfamily.org">Eigen</a> vector.
    vector_type column_vector(int _column,int _tab){
        return vector_type(&(derived()(0,_column,_tab)),this->height());

    }

    const_vector_type column_vector(int _column,int _tab) const{
        return const_vector_type(&(derived()(0,_column,_tab)),this->height());

    }

protected:




    template <class otherDerived,class Op>
    void merge(const DenseLatticeBase<otherDerived> &b,Op op)
    {


        typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;
        typename otherDerived::const_storage_iterator other_it=b.begin();
        for(storage_iterator i=begin(); i<end(); i++)
            *i=to_mdata_type::convert(op(*i,*other_it++));

    }
private:



};



template <class Derived>
template<class otherDerived>
inline
typename DenseProductReturnType<Derived,otherDerived>::type
DenseLatticeBase<Derived>::operator*(const DenseLatticeBase<otherDerived> &b) const
{


    check_mult_dims(b);



    typedef typename DenseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename RType::matrix_type return_matrix_type;
    typedef typename otherDerived::const_matrix_type b_matrix_type;
    RType c(this->height(),b.width(),this->depth());

    for (int i=0; i<this->depth(); i++)
    {


        c.derived().tab_matrix(i)=(tab_matrix(i))*(b.derived().tab_matrix(i));

    }

    return c;
}


template <class Derived>
template<class otherDerived>
inline
typename DenseSolveReturnType<Derived,otherDerived>::type
DenseLatticeBase<Derived>::solve(const DenseLatticeBase<otherDerived> &b) const{


    check_solve_dims(b);



    typedef typename DenseSolveReturnType<Derived,otherDerived>::type RType;

    RType c(this->width(),b.width(),this->depth());

    if (this->width()==this->height()){
        inverse_solve(b,c);
    }
    else if (this->width()<this->height()){
        lsqr_solve(b,c);

    }



    return c;


}


template <class Derived>
template<class otherDerived>
inline void DenseLatticeBase<Derived>::inverse_solve(const DenseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const{

    typedef typename DenseSolveReturnType<Derived,otherDerived>::type RType;
    typedef typename RType::matrix_type return_matrix_type;
    typedef typename otherDerived::const_matrix_type b_matrix_type;


    for (int i=0; i<this->depth(); i++)
    {

        c.derived().tab_matrix(i)=tab_matrix(i).colPivHouseholderQr().solve(b.derived().tab_matrix(i));

    }


}

template <class Derived>
template<class otherDerived>
inline void DenseLatticeBase<Derived>::lsqr_solve(const DenseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const{

    typedef typename DenseSolveReturnType<Derived,otherDerived>::type RType;
    typedef typename RType::matrix_type return_matrix_type;
    typedef typename otherDerived::const_matrix_type b_matrix_type;

    for (int i=0; i<this->depth(); i++)
    {

        c.derived().tab_matrix(i)=tab_matrix(i).jacobiSvd().solve(b.derived().tab_matrix(i));
    }


}


template<class Derived>
void DenseLatticeBase<Derived>::print() const
{
    for (int i=0; i<this->depth(); i++)
    {


        std::cout<< "Tab: " << i << "\n";
        std::cout<< tab_matrix(i)<< "\n";


    }

}

template<class Derived>
void DenseLatticeBase<Derived>::eye()
{
    for (int i=0; i<this->depth(); i++)
    {


        tab_matrix(i).setIdentity();


    }

}

template<class Derived>
template<class otherDerived>
bool DenseLatticeBase<Derived>::operator==(const DenseLatticeBase<otherDerived> &b) const{

    if(this->height()!=b.height()||this->width()!=b.width()||this->depth()!=b.depth())
        return false;

    for(auto it=derived().data_begin(), itb=b.data_begin();it<derived().data_end();++it,++itb)
        if(*it!=*itb)
            return false;

    return true;


}

template<class Derived>
bool DenseLatticeBase<Derived>::save(const std::string & _filename) const
{
    std::ofstream outfile (_filename.c_str(),std::ios_base::binary);
    if (outfile.fail())
        return false;
    uint64_t temp;
    temp=pack754_64((pack754_type)this->height());
    outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));
    temp=pack754_64((pack754_type)this->width());
    outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));
    temp=pack754_64((pack754_type)this->depth());
    outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));

    for(auto it=derived().data_begin();it<derived().data_end();++it){
        temp=pack754_64(*it);
        outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));
    }

    outfile.close();

    return true;

}




/*! @} */

}



#endif // DENSELATTICEBASE_H
