// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef MIA_H
#define MIA_H

#include <array>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


#include "Index.h"
#include "IndexUtil.h"
#include "MIA_Expr.h"
#include "Util.h"

namespace LibMIA
{

namespace internal
{
template<class Derived>
struct data_type<MIA<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<MIA<Derived> >: public index_type<Derived> {};

template<class Derived>
struct order<MIA<Derived> >: public order<Derived> {};

template<class Derived>
struct data_iterator<MIA<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct storage_iterator<MIA<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<MIA<Derived> >: public const_storage_iterator<Derived> {};
}


/** \addtogroup mia Multi-Index Array Classes
*  @{
*/

//!  Base class for multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  multi-index array types. Provides common operations and functions.

  \tparam Derived   should only be DenseMIABase or SparseMIABase type.
*/
template
<

    class Derived
>
class MIA
{

public:

    typedef typename internal::index_type<Derived>::type index_type;
    typedef typename internal::data_type<Derived>::type data_type;
    constexpr static size_t mOrder=internal::order<Derived>::value;
    Derived& derived() { return *static_cast<Derived*>(this); }
    /** \returns a const reference to the derived object */
    const Derived& derived() const { return *static_cast<const Derived*>(this); }

    template<class...Ts>
    auto operator()(Ts...ts)->MIA_Atom<Derived,typename internal::Indicial_Sequence<Ts...>::sequence> {
        static_assert(internal::check_mia_indexing<MIA,Ts...>::type::value,"Number of dimensions must be same as <order> and each given index must be created using the MIAIndex macro.");
        return MIA_Atom<Derived,typename internal::Indicial_Sequence<Ts...>::sequence>(&derived());

    }

    MIA(){
        for(auto & i:m_dims)
            i=0;
        m_dimensionality=0;
    }

    MIA(const std::array<index_type,mOrder > &_dims): m_dims(_dims),m_dimensionality(compute_dimensionality()) {}

    template<typename... Dims>
    MIA(Dims... dims):m_dims{{dims...}},m_dimensionality(compute_dimensionality()) {
        static_assert(internal::check_mia_constructor<MIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
    }

    void init(const std::array<index_type,mOrder> _dims){
        m_dims(_dims);
        m_dimensionality(compute_dimensionality());
    }

    template<class otherDerived>
    void assign(const MIA<otherDerived>& otherMIA,const std::array<typename otherDerived::index_type,internal::order<MIA>::value>& index_order){
        derived().assign(otherMIA,index_order);
    }

//    toLatticeExpression(std::array<size_t> outer_product_indices, std::array<size_t> inner_product_indices,,std::array<size_t> inter_product_indices){
//        return toLatticeCopy(outer_product_indices, inner_product_indices, inter_product_indices );
//
//    }
//
    /** Sets all mia data to one.*/
    void ones(){
        std::fill ( derived().data_begin(), derived().data_end(), 1);
    }

    /** Sets all mia data to one.*/
    void zeros(){
        std::fill ( derived().data_begin(), derived().data_end(), 0);
    }

    /** Sets all mia data to given value.*/
    void init(data_type val){
        std::fill ( derived().data_begin(), derived().data_end(), val);
    }

    index_type dim(size_t i) const{
        assert(i<mOrder);
        return m_dims[i];
    }

    //!  Sets MIA data to uniformly distributed random values.
    /*!
    If MIA is sparse, only elements already designated as nonzero are set to random values.
    Range is specified using parameters. Will throw a MIAParameterException exception if \f$low>high\f$.

    */
    void randu(int low=0, int high=1){
        using namespace boost::numeric;
        if (low>=high){
            throw MIAParameterException("Lower bound of random numbers must be stricly smaller than upper bound.");
        }
        boost::uniform_real<> uni_dist(low,high);
        boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(gen, uni_dist);
        typedef converter<index_type,boost::uniform_real<>::result_type,conversion_traits<index_type,boost::uniform_real<>::result_type>,def_overflow_handler,RoundEven<boost::uniform_real<>::result_type>> to_mdata_type;
        for (auto i=derived().data_begin();i<derived().data_end();++i){
            *i=to_mdata_type::convert(uni());
        }

    }

    const std::array<index_type,internal::order<MIA>::value>& dims() const{
        return m_dims;
    }
    index_type dimensionality() const{
        return m_dimensionality;
    }


    index_type sub2ind(const std::array<index_type,mOrder> & indices){
        return internal::sub2ind(indices,this->dims());
    }

    std::array<index_type,mOrder> ind2sub(index_type idx){
        return internal::ind2sub(idx, this->dims());
    }

    //! Returns scalar data at given indices
    /*!
        \param[in] indices variadic parameter. Will assert a compile error if size of indices!=mOrder or if Indices datatype are not convertible to index_type
    */
    template<typename... Indices>
    const data_type& at(Indices... indices) const {
        static_assert(internal::check_mia_constructor<MIA,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::array<index_type,internal::order<MIA>::value> temp = {{indices...}};
        return at(temp);
    }

    //! Returns scalar data at given indices
    /*!
        \param[in] indices variadic parameter. Will assert a compile error if size of indices!=mOrder or if Indices datatype are not convertible to index_type
    */
    template<typename... Indices>
    data_type& at(Indices... indices) {
        static_assert(internal::check_mia_constructor<MIA,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::array<index_type,internal::order<MIA>::value> temp = {{indices...}};
        return at(temp);
    }


    //! Returns scalar data at given indices
    const data_type& at(const std::array<index_type, mOrder>& indices) const{

        return derived().atIdx(this->sub2ind(indices));

    }

    //! Returns scalar data at given indices
    data_type& at(const std::array<index_type, mOrder> & indices) {

        return derived().atIdx(this->sub2ind(indices));
    }

    //! Performs destructive add (+=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is added like: this->at(x,y,z)+=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,typename index_param_type>
    MIA & plus_equal(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){
        std::plus<data_type> op;
        derived().merge(b.derived(),op,index_order);
        return *this;
    }


    //! Performs destructive subtract (-=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is subtracted like: this->at(x,y,z)-=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,typename index_param_type>
    MIA & minus_equal(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){
        std::minus<data_type> op;
        derived().merge(b.derived(),op,index_order);
        return *this;
    }


protected:

    index_type compute_dimensionality(){
        index_type running_product=1;
        for(auto i=this->m_dims.begin();i<this->m_dims.end();i++){

            running_product*=*i;

        }
        return running_product;
    };

    void init(){
        m_dimensionality=compute_dimensionality();

    }


    std::array<index_type,mOrder> m_dims;
    index_type m_dimensionality;

    template<class otherDerived,class index_param_type>
    void check_merge_dims(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order)
    {
        for(size_t i=0;i<index_order.size();++i)
            if(b.dim(index_order[i])!=m_dims[i])
                throw MIAParameterException("MIA dimensions must be identical for merger operation (+,-, etc).");

    }


};



/*! @} */


} //namespace LibMIA


#endif // MIA_H
