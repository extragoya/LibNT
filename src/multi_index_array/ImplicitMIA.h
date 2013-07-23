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



#include "LibMiaException.h"
#include "Util.h"
#include "IndexUtil.h"
#include "DenseMIABase.h"



//\defgroup
namespace LibMIA
{

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */
namespace internal
{

template<typename T,size_t _order>
struct data_type<ImplicitMIA<T,_order> >
{
    typedef T type;
};

//implicit MIAs dont have references to their data elements
template<typename T,size_t _order>
struct data_type_ref<ImplicitMIA<T,_order> >
{
    typedef T type;
};

template<typename T,size_t _order>
struct const_data_type_ref<ImplicitMIA<T,_order> >
{
    typedef const T type;
};

template<typename T,size_t _order>
struct index_type<ImplicitMIA<T,_order> >
{
    typedef long long type;
};

template<typename T,size_t _order>
struct order<ImplicitMIA<T,_order> >
{
    constexpr static size_t value=_order;
};

//should never be used
template<typename T,size_t _order>
struct data_iterator<ImplicitMIA<T,_order> >
{
    typedef void type;
};

//should never be used
template<typename T,size_t _order>
struct const_data_iterator<ImplicitMIA<T,_order> >
{
    typedef void type;
};

template<typename T,size_t _order>
struct function_type<ImplicitMIA<T,_order> >
{
    typedef typename data_type<ImplicitMIA<T,_order>>::type data_type;
    typedef typename index_type<ImplicitMIA<T,_order>>::type index_type;
    typedef std::function<data_type(index_type)> type;
};



//template<typename T,size_t _order>
//struct const_data_iterator<ImplicitMIA<T,_order> >
//{
//    typedef typename data_type<ImplicitMIA<T,_order>::type data_type;
//    typedef typename index_type<ImplicitMIA<T,_order>::type index_type;
//    typedef typename function_type<ImplicitMIA<T,_order>::type function_type;
//    typedef boost::function_input_iterator<function_type, index_type> type;
//};

template<typename T,size_t _order>
struct FinalDerived<ImplicitMIA<T,_order> >
{
    typedef ImplicitMIA<T,_order> type;
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
*/
template <class _data_type, size_t _order>
class ImplicitMIA: public DenseMIABase<ImplicitMIA<_data_type,_order> >
{





public:

    //! raw data_type
    typedef typename internal::data_type<ImplicitMIA>::type data_type;
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
    ImplicitMIA():MIA<ImplicitMIA<_data_type,_order> >()
    {
        mFunction=this->zero_function();
    }

    //!  Constructs ImplicitMIA of specified size with a given raw data pointer.
    /*!
        \tparam[in] _dims The dimensions size of data. Will assert a compile failure is size is different than _order


    */
    template<class array_index_type>
    ImplicitMIA(const std::array<array_index_type,_order> &_dims):DenseMIABase<ImplicitMIA<_data_type,_order> >(_dims)
    {


        mFunction=this->zero_function();

    }

    //!  Constructs empty ImplicitMIA
    ImplicitMIA(const ImplicitMIA& otherMIA):DenseMIABase<ImplicitMIA<_data_type,_order> >(otherMIA.dims())
    {
        mFunction=otherMIA.mFunction;
    }


    //!  Constructs ImplicitMIA of specified size with a given raw data pointer.
    /*!
        \tparam[in] _dims The dimensions size of data. Will assert a compile failure is size is different than _order


    */
    template<class array_index_type>
    ImplicitMIA(function_type _function,const std::array<array_index_type,_order> &_dims):DenseMIABase<ImplicitMIA<_data_type,_order> >(_dims)
    {


        mFunction=_function;
    }








    //!  Constructs DenseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<typename... Dims>
    ImplicitMIA(Dims... dims):DenseMIABase<ImplicitMIA<_data_type,_order> > {dims...}
    {

        static_assert(internal::check_mia_constructor<ImplicitMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        mFunction=this->zero_function();

    }




    //! Returns size of raw data. For dense cases, this is the same as dimensionality
    std::size_t size() const
    {

        return this->dimensionality();

    }

    //! Returns scalar data at given linear index
    const data_type atIdx(index_type idx) const{

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
    ->decltype(internal::latticeCopy(*this,row_indices,column_indices,tab_indices))
    {
        return this->toLatticeCopy(row_indices,column_indices,tab_indices);

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

    template<class other_data_type>
    DenseMIA<other_data_type,mOrder> make_explicit() const{
        DenseMIA<other_data_type,mOrder> temp(this->dims());
        //get the explicit values
        for(size_t idx=0;idx<this->dimensionality();++idx){

            temp.atIdx(idx)=temp.convert(this->atIdx(idx));
        }
        return temp;

    }

    template<class other_data_type,class index_param_type>
    DenseMIA<other_data_type,mOrder> make_explicit(const std::array<index_param_type,_order>& index_order) const{
        auto new_dims=this->dims();
        internal::reorder_from(this->dims(),index_order,new_dims);
        DenseMIA<other_data_type,mOrder> temp(new_dims);


        index_type curIdx=0;


        for(auto temp_it=temp.data_begin(); temp_it<temp.data_end(); ++temp_it)
        {
            *temp_it=temp.convert(this->atIdx(internal::sub2ind(temp.ind2sub(curIdx++),index_order,this->dims())));

        }
        return temp;

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

private:






};







/*! @} */

}

#endif // IMPLICITMIA_H_INCLUDED
