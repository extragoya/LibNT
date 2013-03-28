// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef SPARSEMIA_H_INCLUDED
#define SPARSEMIA_H_INCLUDED

#include <tuple>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/tuple/tuple.hpp>
//#include <boost/timer/timer.hpp>

#include "tupleit.hh"

#include "LibMiaException.h"
#include "Util.h"
#include "IndexUtil.h"
#include "SparseMIABase.h"


//\defgroup
namespace LibMIA
{

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */
namespace internal
{

template<typename T,size_t _order>
struct data_type<SparseMIA<T,_order> >
{
    typedef T type;
};

template<typename T,size_t _order>
struct index_type<SparseMIA<T,_order> >
{
    typedef long long type;
};

template<typename T,size_t _order>
struct order<SparseMIA<T,_order> >
{
    constexpr static size_t value=_order;
};


template<typename T,size_t _order>
struct Data<SparseMIA<T,_order> >
{
     typedef std::vector<T> type;
};

template<typename T,size_t _order>
struct Indices<SparseMIA<T,_order> >
{
    typedef std::vector<typename index_type<SparseMIA<T,_order> >::type > type;
};

template<typename T,size_t _order>
struct data_iterator<SparseMIA<T,_order> >
{
    typedef typename Data<SparseMIA<T,_order> >::type::iterator type;
};

template<typename T,size_t _order>
struct index_iterator<SparseMIA<T,_order> >
{
    typedef typename Indices<SparseMIA<T,_order> >::type::iterator type;
};

template<typename T,size_t _order>
struct const_index_iterator<SparseMIA<T,_order> >
{
    typedef typename Indices<SparseMIA<T,_order> >::type::const_iterator type;
};

template<typename T,size_t _order>
struct const_data_iterator<SparseMIA<T,_order> >
{
    typedef typename Data<SparseMIA<T,_order> >::type::const_iterator type;
};

template<typename T,size_t _order>
struct full_iterator_tuple<SparseMIA<T,_order> >
{
    typedef typename std::pair<typename data_iterator<SparseMIA<T,_order> >::type,typename index_iterator<SparseMIA<T,_order> >::type> type;
};

template<typename T,size_t _order>
struct const_full_iterator_tuple<SparseMIA<T,_order>>
{
    typedef typename std::pair<typename const_data_iterator<SparseMIA<T,_order> >::type,typename const_index_iterator<SparseMIA<T,_order> >::type> type;

};

template<typename T,size_t _order>
struct full_tuple<SparseMIA<T,_order> >
{
    typedef std::pair<T&,typename index_type<SparseMIA<T,_order>>::type &> type;
};

template<typename T,size_t _order>
struct const_full_tuple<SparseMIA<T,_order> >
{
    typedef std::pair< const T&,const typename index_type<SparseMIA<T,_order> >::type &> type;
};

template<typename T,size_t _order>
struct storage_iterator<SparseMIA<T,_order> >
{
    typedef typename iterators::TupleIt<typename full_iterator_tuple<SparseMIA<T,_order> >::type > type;
};

template<typename T,size_t _order>
struct const_storage_iterator<SparseMIA<T,_order> >
{
    typedef typename iterators::TupleIt<typename const_full_iterator_tuple<SparseMIA<T,_order> >::type > type;
};


}




//!  MIA class for sparse data.
/*!
  Supports addition, multiplication, and solution of, possibly over-determined, systems of
  linear equations. SparseMIA will own underlying raw data.

  \tparam T   the datatype of individual elements.
  \tparam _order   the order (number of indices) of the MIA.
*/
template <class T, size_t _order>
class SparseMIA: public SparseMIABase<SparseMIA<T,_order> >
{





public:

    //! raw data_type
    typedef typename internal::data_type<SparseMIA>::type data_type;
    //! raw index_type
    typedef typename internal::index_type<SparseMIA>::type index_type;
    //! data container type
    typedef typename internal::Data<SparseMIA>::type Data;
    //! index container type
    typedef typename internal::Indices<SparseMIA>::type Indices;
    typedef typename internal::storage_iterator<SparseMIA>::type storage_iterator;
    typedef typename internal::const_storage_iterator<SparseMIA>::type const_storage_iterator;
    //! iterator type for iterating directly through raw data
    typedef typename internal::data_iterator<SparseMIA>::type data_iterator;
    typedef typename internal::index_iterator<SparseMIA>::type index_iterator;
    typedef typename internal::const_data_iterator<SparseMIA>::type const_data_iterator;
    typedef typename internal::const_index_iterator<SparseMIA>::type const_index_iterator;
    typedef typename internal::full_tuple<SparseMIA>::type full_tuple;
    typedef typename internal::const_full_tuple<SparseMIA>::type const_full_tuple;
    //! raw data pointer type
    typedef T* raw_pointer;
    //! order of the MIA
    constexpr static size_t mOrder=_order;


private:
    friend class SparseMIABase<SparseMIA>;
    Data m_data;
    Indices m_indices;


public:





    //!  Constructs empty SparseMIA
    SparseMIA():SparseMIABase<SparseMIA<T,_order> >(), m_data(), m_indices()
    {
    }


    //! Copy constructor for DenseMIAs
    //
    template<class otherMIA, typename boost::enable_if< internal::is_DenseMIA<otherMIA>,int >::type = 0>
    SparseMIA(const otherMIA& denseMIA):SparseMIABase<SparseMIA<T,_order>>(denseMIA.dims())
    {
        this->copy_other_MIA(denseMIA);

    }

    //! Copy constructor for SparseMIAs
    template<class otherMIA, typename boost::enable_if< internal::is_SparseMIA<otherMIA>,int >::type = 0>
    SparseMIA(const otherMIA& sparseMIA):SparseMIABase<SparseMIA<T,_order>>(sparseMIA.dims())
    {

        this->copy_other_MIA(sparseMIA);

    }







    //!  Constructs SparseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<typename... Dims>
    SparseMIA(Dims... dims):SparseMIABase<SparseMIA<T,_order> > {dims...}, m_data(),m_indices()
    {


    }

    //!  Constructs SparseMIA of specified size with a given scalar data and index containers.
    /*!
        Will swap the contents of the container parameters, meaning passed in containers will now be empty, invalidating all previous references, etc.

        \param[in] scalar_data scalar data values
        \param[in] indice_data linear index values - must be same length as scalar_data
        \param[in] _dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than mOrder or if datatype of variables making
                            up _dims are not convertible to index_type

    */
    template<typename... Dims>
    SparseMIA(Data & scalar_data,Indices & indice_data,Dims... dims):SparseMIABase<SparseMIA<T,_order> > {dims...},m_data(),m_indices()
    {

        m_data.swap(scalar_data);
        m_indices.swap(indice_data);

    }

    //!  Constructs SparseMIA of specified size with a given raw scalar data and index data.
    /*!
        Will copy the contents of the raw scalar and index data

        \param[in] scalar_data pointer to scalar data values
        \param[in] indice_data pointer to linear index values - must be same length as scalar_data
        \param[in] data_length number of non-zero elements - ie size of scalar_data and indice_data
        \param[in] _dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than mOrder or if datatype of variables making
                            up _dims are not convertible to index_type

    */
    template<typename... Dims>
    SparseMIA(const data_type * scalar_data,const index_type * & indice_data,size_t _nnz,Dims... dims):SparseMIABase<SparseMIA<T,_order> > {dims...},m_data(scalar_data,scalar_data+_nnz),m_indices(indice_data,indice_data+_nnz)
    {

    }

    //!  Assignment based on given order.
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up

        \param[in] otherMIA the other MIA
        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {2,0,1} this->at(x,y,z)==otherMIA.at(y,z,x).

    */
    template<typename otherDerived,typename index_param_type>
    void assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);


    template<typename otherMIAType,typename boost::enable_if< internal::is_SparseMIA<otherMIAType>,int >::type = 0>
    SparseMIA & operator=(const otherMIAType& otherMIA){
        this->m_dims=otherMIA.dims();
        this->mSortOrder=otherMIA.sort_order();
        this->copy_other_MIA(otherMIA);
        return *this;
    }

    template<typename otherMIAType,typename boost::enable_if< internal::is_DenseMIA<otherMIAType>,int >::type = 0>
    SparseMIA & operator=(const otherMIAType& otherMIA){
        this->m_dims=otherMIA.dims();
        this->copy_other_MIA(otherMIA);
        return *this;
    }


    //! Returns the data container
    const Data & data() const{
        return m_data;
    }
    const Indices& indices() const{
        return m_indices;
    }

    //! Returns a raw pointer to the scalar data
    T* raw_data_ptr(){
        return &m_data[0];
    }

        //! Returns a raw pointer to the scalar data
    const T* raw_data_ptr() const{
        return &m_data[0];
    }

    //! Iterator to the beginning of the raw data
    data_iterator data_begin()
    {
        return m_data.begin();

    }

    //! Iterator to the end of the raw data
    data_iterator data_end()
    {
        return m_data.end();

    }

    //! Iterator to the beginning of the index data
    index_iterator index_begin()
    {
        return m_indices.begin();

    }

    //! Iterator to the end of the index data
    index_iterator index_end()
    {
        return m_indices.end();

    }

    //! Iterator to the beginning of the raw data
    const_data_iterator data_begin() const
    {
        return m_data.begin();

    }

    //! Iterator to the end of the raw data
    const_data_iterator data_end() const
    {
        return m_data.end();

    }

    //! Iterator to the beginning of the index data
    const_index_iterator index_begin() const
    {
        return m_indices.begin();

    }

    //! Iterator to the end of the index data
    const_index_iterator index_end() const
    {
        return m_indices.end();

    }



    //! Converts a scalar value to data_type
    /*!
        \tparam from_data_type the data_type you are converting from
    */
    template<class from_data_type>
    data_type convert(const from_data_type from) const{
        using namespace boost::numeric;
        typedef converter<data_type,from_data_type,conversion_traits<data_type,from_data_type>,def_overflow_handler,RoundEven<from_data_type>> to_mdata_type;
        return to_mdata_type::convert(from);
    }

    //! Returns size of raw data. For sparse cases, this is the number of nonzeros
    std::size_t size() const
    {

        return m_data.size();

    }

    void clear()
    {

        m_data.clear();
        m_indices.clear();

    }



protected:

    template<typename otherDerived, typename Op,typename index_param_type>
    void scanMerge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order);

    template<typename otherDerived, typename Op,typename index_param_type>
    void sortMerge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order);
private:






};

template <class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type>
void  SparseMIA<T,_order>::scanMerge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order)
{

    //get the order of lhs indices in terms of rhs
    auto lhsOrder=internal::reverseOrder(index_order);
    //we also need to reorder b's sort order (which may not be {0,1,2, etc.}) using the index order
    lhsOrder= internal::reOrderArray(b.sort_order(), lhsOrder);
    //print_array(lhsOrder,"lhsOrder");
    //boost::timer::cpu_timer sort_t;

    //this->change_sort_order(lhsOrder);

    //std::sort(this->index_begin(),this->index_end());
    //std::sort(this->data_begin(),this->data_end());
    if(lhsOrder!=this->mSortOrder)
       this->sort(lhsOrder);
    //std::cout << "Scan sort " << boost::timer::format(sort_t.elapsed()) << std::endl;
    //this->print();
    size_t old_size=this->size();
    this->resize(this->size()+b.size());
    //boost::timer::cpu_timer merge_t;
    auto new_end=internal::merge_sparse_storage_containers(this->storage_begin(),this->storage_begin()+old_size,b.storage_begin(),b.storage_end(),op);
    //std::cout << "Scan merge " << boost::timer::format(merge_t.elapsed()) << std::endl;
    this->resize(new_end-this->storage_begin());

}

template <class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type>
void  SparseMIA<T,_order>::sortMerge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order)
{

    //get the order of lhs indices in terms of rhs
    auto lhsOrder=internal::reverseOrder(index_order);
    //we also need to reorder b's sort order (which may not be {0,1,2, etc.}) using the index order
    lhsOrder= internal::reOrderArray(b.sort_order(), lhsOrder);
    if(lhsOrder!=this->mSortOrder)
       this->change_sort_order(lhsOrder);
    //this->print();
    size_t old_size=this->size();
    this->resize(this->size()+b.size());
    this->mIsSorted=false;
    if(boost::is_same<std::minus<data_type>,Op>::value){
        //since stable sort is slower (inherently and also because we have to use tupleit with stable_sort), we just negate b's datatype
        //during the copy process, and then just perform addition
        std::function<data_type(typename SparseMIABase<otherDerived>::data_type)> copy_function = [&](const typename SparseMIABase<otherDerived>::data_type & _other_data)
        {
            return this->convert(-1*_other_data);
        };
        std::copy(boost::make_transform_iterator(b.data_begin(), copy_function),
                  boost::make_transform_iterator(b.data_end(), copy_function),
                  this->data_begin()+old_size);
        std::copy(b.index_begin(),b.index_end(),this->index_begin()+old_size);
        this->collect_duplicates(std::plus<data_type>());

    }
    else{
        std::copy(b.storage_begin(),b.storage_end(),this->storage_begin()+old_size);
        this->collect_duplicates(op);
    }




}



template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void SparseMIA<T,_order>::assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{

    std::cout <<"We got to sparse assign" << std::endl;
    std::array<size_t,mOrder> converted_index_order=array_converter<size_t>::convert(index_order);
    internal::reorder_from(otherMIA.dims(),index_order,this->m_dims);
    internal::reorder_to(otherMIA.sort_order(),index_order,this->mSortOrder);
    if(otherMIA.sort_order()!=this->mSortOrder)
        this->mIsSorted=false;
    else
        this->mIsSorted=otherMIA.is_sorted();

    this->resize(otherMIA.size());

    auto otherIt=otherMIA.storage_begin();
    std::for_each(this->storage_begin(),this->storage_end(),[&](full_tuple & cur_tuple){
        cur_tuple=this->convert(*(otherIt++));
    });




}
//
//template<class T, size_t _order>
//template<typename otherDerived>
//SparseMIA<T,_order>& SparseMIA<T,_order>::operator=(const SparseMIABase<otherDerived>& otherMIA)
//{
//
//    std::cout <<"We got to sparse =" << std::endl;
//    this->m_dims=otherMIA.dims();
//    this->mSortOrder=otherMIA.sort_order();
//    this->mIsSorted=otherMIA.is_sorted();
//    otherMIA.print();
//    this->resize(otherMIA.size());
//
//    auto otherIt=otherMIA.storage_begin();
//    std::for_each(this->storage_begin(),this->storage_end(),[this,&otherIt](full_tuple cur_tuple){
//        this->data_val(cur_tuple)=this->convert(boost::get<0>(*otherIt));
//        this->index_val(cur_tuple)=boost::get<1>(*otherIt++);
//    });
//
//    return *this;
//
//
//
//}


/*! @} */

}

#endif // SparseMIA_H_INCLUDED
