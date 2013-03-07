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
#include <boost/utility/enable_if.hpp>

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
    typedef typename std::tuple<typename data_iterator<SparseMIA<T,_order> >::type,typename index_iterator<SparseMIA<T,_order> >::type> type;
};

template<typename T,size_t _order>
struct const_full_iterator_tuple<SparseMIA<T,_order>>
{
    typedef typename std::tuple<typename const_data_iterator<SparseMIA<T,_order> >::type,typename const_index_iterator<SparseMIA<T,_order> >::type> type;

};

template<typename T,size_t _order>
struct full_tuple<SparseMIA<T,_order> >
{
    typedef std::tuple<T&,typename index_type<SparseMIA<T,_order>>::type &> type;
};

template<typename T,size_t _order>
struct const_full_tuple<SparseMIA<T,_order> >
{
    typedef std::tuple< const T&,const typename index_type<SparseMIA<T,_order> >::type &> type;
};

template<typename T,size_t _order>
struct storage_iterator<SparseMIA<T,_order> >
{
    typedef typename iterators::TupleIt<typename full_iterator_tuple<SparseMIA<T,_order> >::type > type;
};

template<typename T,size_t _order>
struct const_storage_iterator<SparseMIA<T,_order> >
{
    typedef typename boost::zip_iterator<typename const_full_iterator_tuple<SparseMIA<T,_order> >::type > type;
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
    //! raw data pointer type
    typedef T* raw_pointer;
    //! order of the MIA
    constexpr static size_t mOrder=_order;


private:

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
        //count the number of nnzs
        size_t nnz=0;
        for(auto it=denseMIA.data_begin();it<denseMIA.data_end();++it)
            if(*it)
                ++nnz;

        //allocate the required data
        m_data.resize(nnz);
        m_indices.resize(nnz);

        //set m_data and m_indices
        index_type idx=0;
        nnz=0;
        for(auto it=denseMIA.data_begin();it<denseMIA.data_end();++it){
            if(*it){
                m_data[nnz]=*it;
                m_indices[nnz++]=idx;
            }
            ++idx;
        }

    }


//    //!  Copy constructor.
//    /*!
//
//
//    */
//    SparseMIA(const SparseMIA& otherMIA):SparseMIABase<SparseMIA<T,_order> >(otherMIA.dims())
//    {
//
//
//        m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
//        m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
//        std::copy(otherMIA.data_begin(),otherMIA.data_end(),this->data_begin());
//
//
//    }
//
//    //!  Copy constructor.
//    /*!
//        If otherMIA's datatype is different than this->data_type, then individual entries will be converted.
//
//    */
//    template<class otherDerived>
//    SparseMIA(const SparseMIABase<otherDerived>& otherMIA):SparseMIABase<SparseMIA<T,_order> >(otherMIA.dims())
//    {
//
//
//
//        static_assert(internal::order<otherDerived>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");
//        m_smart_raw_ptr.reset(new T[this->m_dimensionality]);
//        m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
//        size_t idx=0;
//        std::for_each( otherMIA.data_begin(), otherMIA.data_end(), [this,&idx] (typename otherDerived::data_type val)
//        {
//            *(this->data_begin()+idx++)=this->convert<typename otherDerived::data_type>(val);
//        } );
//
//
//
//
//    }


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

//    //!  Assignment based on given order.
//    /*!
//
//        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
//        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up
//
//        \param[in] otherMIA the other MIA
//        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {3,1,2} this->at(1,2,3)==otherMIA.at(2,3,1).
//                                Will assert a compile failure is size of index_order is not the same as this->mOrder
//    */
//    template<typename otherDerived,typename index_param_type>
//    void assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);
//
//    //!  Assignment operator.
//    /*!
//        \param[in] otherMIA
//        If the data_type of otherMIA is not the same as this, the scalar data will be converted.
//    */
//    template<typename otherDerived>
//    SparseMIA& operator=(const SparseMIABase<otherDerived>& otherMIA);
//
//    //!  Straight-out assignment.
//    SparseMIA& operator=(const SparseMIA& otherMIA){
//        return *this=static_cast<const SparseMIABase<SparseMIA>&>(otherMIA);
//    }

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
        return boost::numeric::converter<data_type,from_data_type>::convert(from);
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

private:






};


//template<class T, size_t _order>
//template<typename otherDerived>
//SparseMIA<T,_order>& SparseMIA<T,_order>::operator=(const SparseMIABase<otherDerived>& otherMIA)
//{
//    static_assert(_order==internal::order<otherDerived>::value,"Orders of MIAs must be the same to be assigned");
//
//
//    if(this->m_dimensionality!=otherMIA.dimensionality()){
//        if(hasOwnership){
//            this->m_dimensionality=otherMIA.dimensionality();
//            this->m_dims=otherMIA.dims();
//            smart_raw_pointer temp_ptr(new T[this->m_dimensionality]);
//            m_smart_raw_ptr.swap(temp_ptr);
//            temp_ptr.release();
//            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
//        }
//        else
//            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
//    }
//    else if(this->m_dims!=otherMIA.dims())
//        this->m_dims=otherMIA.dims();
//
//    typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;
//    for(auto it1=this->data_begin(),it2=otherMIA.data_begin();it1<this->data_end();++it1,++it2)
//        *it1=to_mdata_type::convert(*it2);
//
//    return *this;
//
//}
//
//
//
//
//
//
//template<class T, size_t _order>
//template<typename otherDerived,typename index_param_type>
//void SparseMIA<T,_order>::assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
//{
//    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");
//
//
//    if(this->m_dimensionality!=otherMIA.dimensionality()){
//        if(hasOwnership){
//            internal::collect_dimensions_from_order(otherMIA,index_order,this->m_dims);
//            this->m_dimensionality=otherMIA.dimensionality();
//            smart_raw_pointer temp_ptr(new T[this->m_dimensionality]);
//            m_smart_raw_ptr.swap(temp_ptr);
//            temp_ptr.release();
//            m_Data.reset(new data_container_type(m_smart_raw_ptr.get(),this->m_dims,boost::fortran_storage_order()));
//        }
//        else
//            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
//    }
//    else
//        internal::collect_dimensions_from_order(otherMIA,index_order,this->m_dims);
//
//
//    index_type curIdx=0;
//
//    auto other_it=otherMIA.data_begin();
//    for(auto this_it=this->data_begin(); this_it<this->data_end(); ++this_it)
//    {
//        *this_it=this->convert(*(other_it+sub2ind(ind2sub(curIdx++, this->dims()),index_order,otherMIA.dims())));
//
//    }
//
//
//}




/*! @} */

}

#endif // SparseMIA_H_INCLUDED
