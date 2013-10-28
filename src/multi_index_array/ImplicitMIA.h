// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef IMPLICITMIA_H_INCLUDED
#define IMPLICITMIA_H_INCLUDED

#include <type_traits>
#include <iostream>
#include <algorithm>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include "LibMiaException.h"
#include "LibMIAUtil.h"
#include "IndexUtil.h"
#include "DenseMIABase.h"



//\defgroup
namespace LibMIA
{

//forward declaration of data iterator
template<typename ImplictMIAType>
class implicit_iter;



/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */
namespace internal
{

template<typename T,size_t _order,bool isRef>
struct data_type<ImplicitMIA<T,_order, isRef> >
{
    typedef T type;
};

//implicit MIAs dont have references to their data elements
template<typename T,size_t _order>
struct data_type_ref<ImplicitMIA<T,_order, false> >
{
    typedef T type;
};

template<typename T,size_t _order>
struct const_data_type_ref<ImplicitMIA<T,_order, false> >
{
    typedef const T type;
};

//implicit MIAs dont have references to their data elements
template<typename T,size_t _order>
struct data_type_ref<ImplicitMIA<T,_order, true> >
{
    typedef T & type;
};

template<typename T,size_t _order>
struct const_data_type_ref<ImplicitMIA<T,_order, true> >
{
    typedef const T & type;
};

template<typename T,size_t _order,bool isRef>
struct index_type<ImplicitMIA<T,_order, isRef> >
{
    typedef long long type;
};

template<typename T,size_t _order,bool isRef>
struct order<ImplicitMIA<T,_order, isRef> >
{
    constexpr static size_t value=_order;
};

//should never be used
template<typename T,size_t _order,bool isRef>
struct data_iterator<ImplicitMIA<T,_order, isRef> >
{
    typedef implicit_iter<ImplicitMIA<T,_order, isRef>> type;
};

//should never be used
template<typename T,size_t _order,bool isRef>
struct const_data_iterator<ImplicitMIA<T,_order, isRef> >
{
    typedef implicit_iter<const ImplicitMIA<T,_order, isRef>> type;
};

template<typename T,size_t _order,bool isRef>
struct function_type<ImplicitMIA<T,_order, isRef> >
{
    typedef typename data_type_ref<ImplicitMIA<T,_order, isRef>>::type data_type_ref;
    typedef typename index_type<ImplicitMIA<T,_order, isRef>>::type index_type;
    typedef std::function<data_type_ref(const index_type &)> type;
};



//template<typename T,size_t _order>
//struct const_data_iterator<ImplicitMIA<T,_order> >
//{
//    typedef typename data_type<ImplicitMIA<T,_order>::type data_type;
//    typedef typename index_type<ImplicitMIA<T,_order>::type index_type;
//    typedef typename function_type<ImplicitMIA<T,_order>::type function_type;
//    typedef boost::function_input_iterator<function_type, index_type> type;
//};

template<typename T,size_t _order,bool isRef>
struct FinalDerived<ImplicitMIA<T,_order, isRef> >
{
    typedef ImplicitMIA<T,_order, isRef> type;
};

}




//!  MIA class for implicit MIAs (where values are defined by a function). Can be treated as a read-only type of MIA
/*!
  Supports addition, multiplication, and solution of, possibly over-determined, systems of
  linear equations.

  If the data-generating function is provided by the user, care must be taken to ensure that any objects referred to by the function
  (if any) are destroyed prior to the function call.

  \tparam T   the datatype of individual elements.
  \tparam _order   the order (number of indices) of the MIA.
  \tparam isRef if the ImplicitMIA references another Dense MIA raw data, set to true
*/
template <class _data_type, size_t _order,bool isRef>
class ImplicitMIA: public DenseMIABase<ImplicitMIA<_data_type,_order,isRef> >
{





public:

    //! raw data_type
    typedef typename internal::data_type<ImplicitMIA>::type data_type;
    //! raw data_type reference
    typedef typename internal::data_type_ref<ImplicitMIA>::type data_type_ref;
    //! const raw data_type reference
    typedef typename internal::const_data_type_ref<ImplicitMIA>::type const_data_type_ref;

    //! raw index_type
    typedef typename internal::index_type<ImplicitMIA>::type index_type;

//    //! iterator type for iterating directly through raw data
//    typedef typename internal::const_data_iterator<ImplicitMIA>::type const_data_iterator;

    //! implict function type
    typedef typename internal::function_type<ImplicitMIA>::type function_type;

    //! final derived type
    typedef typename internal::FinalDerived<ImplicitMIA>::type FinalDerived;

    //! order of the MIA
    constexpr static size_t mOrder=_order;

    typedef typename internal::data_iterator<ImplicitMIA>::type data_iterator;
    typedef typename internal::const_data_iterator<ImplicitMIA>::type const_data_iterator;


private:

    function_type mFunction;


public:



    FinalDerived& final_derived()
    {
        return *this;
    }
    /** \returns a const reference to the derived object */
    const FinalDerived& final_derived() const
    {
        return *this;
    }

    //!  Constructs empty ImplicitMIA
    ImplicitMIA():DenseMIABase<ImplicitMIA<_data_type,_order,isRef> >()
    {
        mFunction=this->zero_function();
    }

    //!  Constructs ImplicitMIA of specified size with a given raw data pointer.
    /*!
        \tparam[in] _dims The dimensions size of data. Will assert a compile failure is size is different than _order


    */
    template<class array_index_type>
    ImplicitMIA(const std::array<array_index_type,_order> &_dims):DenseMIABase<ImplicitMIA<_data_type,_order,isRef> >(_dims)
    {


        mFunction=this->zero_function();

    }

    //!  Constructs empty ImplicitMIA
    ImplicitMIA(const ImplicitMIA& otherMIA):DenseMIABase<ImplicitMIA<_data_type,_order,isRef> >(otherMIA.dims())
    {
        mFunction=otherMIA.mFunction;
    }


    //!  Constructs ImplicitMIA of specified size with a given raw data pointer.
    /*!
        \tparam[in] _dims The dimensions size of data. Will assert a compile failure is size is different than _order


    */
    template<class array_index_type>
    ImplicitMIA(const function_type &_function,const std::array<array_index_type,_order> &_dims):DenseMIABase<ImplicitMIA<_data_type,_order,isRef> >(_dims)
    {


        mFunction=_function;
    }


    //!  Constructs DenseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<typename... Dims>
    ImplicitMIA(const function_type &_function,Dims... dims):DenseMIABase<ImplicitMIA<_data_type,_order,isRef> > {dims...}
    {

        static_assert(internal::check_mia_constructor<ImplicitMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        mFunction=_function;

    }

    ImplicitMIA& operator=(const ImplicitMIA & otherMIA){
        this->m_dims=otherMIA.dims();
        this->mFunction=otherMIA.get_function();
        this->m_dimensionality=otherMIA.dimensionality();
        return *this;
    }



    //!  Constructs DenseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<typename... Dims>
    ImplicitMIA(Dims... dims):DenseMIABase<ImplicitMIA<_data_type,_order,isRef> > {dims...}
    {

        static_assert(internal::check_mia_constructor<ImplicitMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        mFunction=this->zero_function();

    }


    template<typename Derived>
    ImplicitMIA & operator=(const DenseMIABase<Derived> & otherMIA)
    {

        check_if_assign_legal<isRef>();
        if(this->dims()!=otherMIA.dims()){
            throw MIAParameterException("To assign to an ImplicitMIA that refers to another's underlying data, dims of two  MIAs must be identical.");
        }
        for(size_t idx=0;idx<this->dimensionality();++idx){
            this->atIdx(idx)=this->convert(otherMIA.atIdx(idx));
        }
        return *this;
    }




    //! Returns size of raw data. For dense cases, this is the same as dimensionality
    std::size_t size() const
    {

        return this->dimensionality();

    }

    //! Returns scalar data at given linear index
    inline const_data_type_ref atIdx(index_type idx) const{

        //return lin index
        const_data_type_ref ret= mFunction(idx);
        return ret;
    }

    //! Returns scalar data at given linear index
    inline data_type_ref atIdx(index_type idx) {

        //return lin index
        return mFunction(idx);
    }


    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom when the MIA in question is a temp object.
    /*!

    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    auto toLatticeDiscard(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices)
    ->decltype(internal::latticeCopy(*this,row_indices,column_indices,tab_indices))
    {
        return this->toLatticeCopy(row_indices,column_indices,tab_indices);

    }

    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom when the MIA in question is a temp object.
    /*!

    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    auto toLatticeExpression(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices)
    ->decltype(this->toLatticeCopy(row_indices,column_indices,tab_indices)) const
    {
        return this->toLatticeCopy(row_indices,column_indices,tab_indices);

    }


    auto toStraightLattice(size_t number_of_row_indices, size_t number_of_column_indices)
    ->decltype(this->toStraightLatticeCopy(number_of_row_indices,number_of_column_indices)) const
    {

        return this->toStraightLatticeCopy(number_of_row_indices,number_of_column_indices);

    }

//    template<class otherDerived,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
//    typename MIAMergeReturnType<ImplicitMIA,otherDerived>::type  plus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{
//        std::plus<data_type> op;
//        return implicit_merge(b.derived(),op,index_order);
//    }
//
//    template<class otherDerived,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
//    typename MIAMergeReturnType<ImplicitMIA,otherDerived>::type  minus_(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const{
//        std::minus<data_type> op;
//        return implicit_merge(b.derived(),op,index_order);
//    }

    function_type & get_function(){
        return mFunction;
    }

    const function_type & get_function() const{
        return mFunction;
    }

    template<class other_data_type=data_type>
    DenseMIA<other_data_type,mOrder> make_explicit() const{
        DenseMIA<other_data_type,mOrder> temp(this->dims());
        //get the explicit values
        //get faster execution if I make different code structures instead of using opemmp's if statement
        if(this->dimensionality()>=PARALLEL_TOL){
            #pragma omp parallel for
            for(size_t idx=0;idx<this->dimensionality();++idx){

                temp.atIdx(idx)=temp.convert(this->atIdx(idx));
            }
        }
        else
        {
            for(size_t idx=0;idx<this->dimensionality();++idx){

                temp.atIdx(idx)=temp.convert(this->atIdx(idx));
            }
        }
        return temp;

    }



    template<class other_data_type,class index_param_type>
    DenseMIA<other_data_type,mOrder> make_explicit(const std::array<index_param_type,_order>& index_order) const{
        auto temp_storage=this->dims();
        internal::reorder_from(this->dims(),index_order,temp_storage);
        DenseMIA<other_data_type,mOrder> temp(temp_storage);





        auto dim_accumulator=internal::createDimAccumulator(temp_storage,index_order); //precompute the demoninators needed to convert from linIdx to a full index, using new_dims
        //get faster execution if I make different code structures instead of using opemmp's if statement
        if(this->dimensionality()>=PARALLEL_TOL){
            #pragma omp parallel for
            for(auto temp_it=temp.data_begin(); temp_it<temp.data_end(); ++temp_it)
            {

                *temp_it=temp.convert(this->atIdx(internal::getShuffleLinearIndex((index_type)(temp_it-temp.data_begin()),this->dims(),dim_accumulator)));

            }
        }
        else
        {
            for(auto temp_it=temp.data_begin(); temp_it<temp.data_end(); ++temp_it)
            {

                *temp_it=temp.convert(this->atIdx(internal::getShuffleLinearIndex((index_type)(temp_it-temp.data_begin()),this->dims(),dim_accumulator)));

            }
        }

        return temp;

    }

    data_iterator data_begin(){
        return data_iterator(this,0);
    }

    data_iterator data_end(){
        return data_iterator(this,this->dimensionality());
    }

    const_data_iterator data_begin() const{
        return const_data_iterator(this,0);
    }

    const_data_iterator data_end() const{
        return const_data_iterator(this,this->dimensionality());
    }



protected:

    function_type zero_function(){
        function_type zero_func=[](index_type dummy){
            return 0;
        };
        return zero_func;

    }

    template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
    typename MIAMergeReturnType<ImplicitMIA,otherDerived>::type
    implicit_merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,mOrder>& index_order) const
    {





        this->check_merge_dims(b,index_order);
        static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

        return internal::perform_implicit_merge(*this, b,op,index_order);


    }

    //If ImplicitMIAs refer to some DenseMIA data held elsewhere, ie isRef is true, then under strict conditions assignment is legal. Otherwise, assignment is always illegal. These two functions
    //ensure that a compilation error will trigger if someone tries to assign with isRef set to false.
    template<bool should_do_assign,typename boost::disable_if_c<should_do_assign,int>::type=0>
    void check_if_assign_legal()
    {
        //use delayed parsing, so the static assert will only trigger if the function is actually used within a compilation unit
        struct fake : std::false_type{};
        static_assert(fake::value,"ImplictMIA must have its isRef tparam set to true to perform assignment");
    }

    template<bool should_do_assign,typename boost::enable_if_c<should_do_assign,int>::type=0>
    void check_if_assign_legal()
    {

    }

private:






};




//!Iterator for ImplicitMIAs.
/*!
        When ImplicitMIA's isRef is false, just returns data_type upon dereference, otherwise returns &data_type. In the latter case, interoperability
        with std::algorithms has not been tested for swapping and copying algorithms, like std::sort. Can handle const, by passing a const ImplicitMIA
        type as ImplicitMIAType (typedefs in the traits class)
*/
template<class ImplicitMIAType>
class implicit_iter
  : public boost::iterator_facade<
       implicit_iter<ImplicitMIAType>,
       typename boost::mpl::if_<
            boost::is_const<ImplicitMIAType>,
            const typename internal::data_type<typename std::remove_const<ImplicitMIAType>::type>::type,
            typename internal::data_type<ImplicitMIAType>::type
        >::type,
        boost::random_access_traversal_tag,
        typename boost::mpl::if_<
            boost::is_const<ImplicitMIAType>,
            const typename internal::data_type_ref<typename std::remove_const<ImplicitMIAType>::type>::type, //is isRef is false in ImplicitMIAType, will return data_type instead of data_type&
            typename internal::data_type_ref<ImplicitMIAType>::type
        >::type
    >
{
private:
    typedef ImplicitMIAType miaType;
    typedef typename boost::iterator_difference<implicit_iter>::type difference_type;
    typedef typename boost::iterator_reference<implicit_iter>::type reference;

    friend class boost::iterator_core_access;
    template <class> friend class implicit_iterator;

    typedef typename internal::index_type<miaType>::type index_type;

    miaType * mia_ref; //pointer to the mia and the current idx is how the iterator keeps track of how to return data
    index_type mIdx;

    struct enabler {};

    template <class OtherMIAType>
    bool equal(implicit_iter<OtherMIAType> const& other) const
    {
        return this->mia_ref == other.mia_ref && this->mIdx==other.mIdx;
    }

    void increment()
    { mIdx++; }

    void decrement()
    { mIdx--; }

    void advance(difference_type _diff){
        mIdx+=_diff;
    }

    template <class OtherMIAType>
    difference_type distance_to(implicit_iter<OtherMIAType> const& other) const{
        return other.mIdx-this->mIdx;
    }

    reference dereference() const
    { return mia_ref->atIdx(mIdx); }

public:


    implicit_iter()
      : mia_ref(0),mIdx(0) {}

    implicit_iter(miaType * _mia)
      : mia_ref(_mia),mIdx(0) {}

    explicit implicit_iter(miaType * _mia, index_type _idx)
      : mia_ref(_mia),mIdx(_idx) {}

    template <class OtherMIAType>
    implicit_iter(
      implicit_iter<OtherMIAType> const& other
    , typename boost::enable_if<
          boost::is_convertible<OtherMIAType*,ImplicitMIAType*>
        , enabler
      >::type = enabler()
  )
    : mia_ref(other.mia_ref),mIdx(other.mIdx) {}



};



/*! @} */

}

#endif // IMPLICITMIA_H_INCLUDED
