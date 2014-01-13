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


#include "tupleit.hh"

#include "LibMiaException.h"
#include "LibMIAUtil.h"
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
struct data_type_ref<SparseMIA<T,_order> >
{
    typedef T & type;
};
template<typename T,size_t _order>
struct const_data_type_ref<SparseMIA<T,_order> >
{
    typedef const T & type;
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


template<typename T,size_t _order>
struct FinalDerived<SparseMIA<T,_order> >
{
    typedef SparseMIA<T,_order> type;
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
    typedef typename internal::FinalDerived<SparseMIA>::type FinalDerived;
    //! raw data pointer type
    typedef T* raw_pointer;
    //! order of the MIA
    constexpr static size_t mOrder=_order;


private:
    friend class SparseMIABase<SparseMIA>;
    //!container for data entries - expected to allow Random Access Iterators
    Data m_data;
    //!container for index entries - expected to allow Random Access Iterators
    Indices m_indices;


public:

    FinalDerived& final_derived() {

        return *this;
    }
    /** \returns a const reference to the derived object */
    const FinalDerived& final_derived() const {

        return *this;
    }



    //!  Constructs empty SparseMIA
    SparseMIA():SparseMIABase<SparseMIA<T,_order> >(), m_data(), m_indices()
    {
    }


    //! Copy constructor for DenseMIAs
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

    //! Copy constructor for SparseMIAs
    SparseMIA(const SparseMIA& sparseMIA):SparseMIABase<SparseMIA<T,_order>>(sparseMIA.dims())
    {

        this->copy_other_MIA(sparseMIA);

    }



    SparseMIA(SparseMIA && otherMIA):SparseMIABase<SparseMIA<T,_order>>(otherMIA.dims())
    {


        this->setLinIdxSequence(otherMIA.linIdxSequence());
        this->setSorted(otherMIA.is_sorted());
        m_data.swap(otherMIA.m_data);
        m_indices.swap(otherMIA.m_indices);
        this->mSolveInfo=otherMIA.solveInfo();

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

    //!  Constructs SparseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims array parameter to specify size. Will assert a compile failure if array_index_type of variables making
                            up _dims are not convertible to index_type

    */
    template<class array_index_type>
    SparseMIA(const std::array<array_index_type,_order> &_dims):SparseMIABase<SparseMIA<T,_order> > (_dims), m_data(),m_indices()
    {


    }


    //!  Constructs SparseMIA of specified size with a given scalar data and index containers.
    /*!
        Will swap the contents of the container parameters, meaning passed in containers will now be empty, invalidating all previous references, etc.
        If size of data and index containers doesn't equal the dimensionality of the passed-in dimensions, will throw an MIAParameterException.

        \param[in] scalar_data scalar data values
        \param[in] indice_data linear index values - must be same length as scalar_data
        \param[in] _dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than mOrder or if datatype of variables making
                            up _dims are not convertible to index_type

    */
    template<typename... Dims>
    SparseMIA(Data && scalar_data,Indices && indice_data,Dims... dims):SparseMIABase<SparseMIA<T,_order> > {dims...},m_data(),m_indices()
    {

        if(scalar_data.size()!=indice_data.size()){
            throw MIAParameterException("SparseMIA Constructor: Data and index container parameters must have equal size.");
        }

        m_data.swap(scalar_data);
        m_indices.swap(indice_data);

    }

    //!  Constructs SparseMIA of specified size with a given scalar data and index containers.
    /*!
        Will swap the contents of the container parameters, meaning passed in containers will now be empty, invalidating all previous references, etc.
        If size of data and index containers doesn't equal the dimensionality of the passed-in dimensions, will throw an MIAParameterException.

        \param[in] scalar_data scalar data values
        \param[in] indice_data linear index values - must be same length as scalar_data
        \param[in] _dims array parameter to specify size. Will assert a compile failure if array_index_type of variables making
                            up _dims are not convertible to index_type

    */
    template<typename array_index_type>
    SparseMIA(Data && scalar_data,Indices && indice_data,const std::array<array_index_type,_order> &_dims,bool _is_sorted=true):SparseMIABase<SparseMIA<T,_order> >(_dims,_is_sorted),m_data(),m_indices()
    {

        if(scalar_data.size()!=indice_data.size()){
            throw MIAParameterException("SparseMIA Constructor: Data and index container parameters must have equal size.");
        }


        m_data.swap(scalar_data);
        m_indices.swap(indice_data);

    }




    //!  Constructs SparseMIA of specified size with a given SparseLattice
    /*!
        Will swap the contents of the SparseLattice, meaning passed in containers will now be empty, invalidating all previous references, etc.

        \param[in] _dims array parameter to specify size. Will assert a compile failure if array_index_type of variables making
                            up _dims are not convertible to index_type
        \param[in] _lat an rvalue reference to a SparseLattice. Do not use the SparseLattice reference after calling this constructor

    */
    template<class array_index_type>
    SparseMIA(const std::array<array_index_type,_order> &_dims,SparseLattice<data_type>&& _lat):SparseMIABase<SparseMIA<T,_order> >(_dims,_lat.is_sorted()&&_lat.linIdxSequence()==ColumnMajor),m_data(),m_indices()
    {


        _lat.change_sort_order(ColumnMajor); //make sure lattice is columnmajor, so that indices are calculated in the expected linIdxSequence
        m_data.swap(_lat.data());
        m_indices.swap(_lat.indices());


    }

    //!  Constructs SparseMIA of specified size with a given raw scalar data and index data.
    /*!
        Will copy the contents of the raw scalar and index data.

        \param[in] scalar_data pointer to scalar data values
        \param[in] indice_data pointer to linear index values - must be same length as scalar_data
        \param[in] data_length number of non-zero elements - ie size of scalar_data and indice_data
        \param[in] _dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than mOrder or if datatype of variables making
                            up _dims are not convertible to index_type

    */
    template<typename... Dims>
    SparseMIA(const data_type * scalar_data,const index_type * indice_data,size_t _nnz,Dims... dims):SparseMIABase<SparseMIA<T,_order> > {dims...},m_data(scalar_data,scalar_data+_nnz),m_indices(indice_data,indice_data+_nnz)
    {

    }




    //!  Assignment based on given order.
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up. Even otherMIA is
        sorted, *this may not be if index_order designates a shuffle.

        \param[in] otherMIA the other MIA
        \param[in] index_order An index reshuffle denoting how otherMIA is indexed based on *this. E.g., if order is {2,0,1} this->at(0,1,2)==otherMIA.at(1,2,0).

    */
    template<typename otherDerived,typename index_param_type>
    void assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);


    //!  Assignment based on given order from a Dense MIA.
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up. Even otherMIA is
        sorted, *this may not be if index_order designates a shuffle.

        \param[in] otherMIA the other MIA
        \param[in] index_order An index reshuffle denoting how otherMIA is indexed based on *this. E.g., if order is {2,0,1} this->at(0,1,2)==otherMIA.at(1,2,0).

    */
    template<typename otherDerived,typename index_param_type>
    void assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);



    //!  Sparse assignment
    template<typename otherMIAType,typename boost::enable_if< internal::is_SparseMIA<otherMIAType>,int >::type = 0>
    SparseMIA & operator=(const otherMIAType& otherMIA){

        this->mLinIdxSequence=otherMIA.linIdxSequence();
        this->copy_other_MIA(otherMIA);
        return *this;
    }
    //!  Dense assignment
    template<typename otherMIAType,typename boost::enable_if< internal::is_DenseMIA<otherMIAType>,int >::type = 0>
    SparseMIA & operator=(const otherMIAType& otherMIA){


        this->copy_other_MIA(otherMIA);
        return *this;
    }
    SparseMIA& operator=(const SparseMIA& otherMIA) {



        this->copy_other_MIA(otherMIA);
        return *this;

    }

    //!  Sparse move assignment
    SparseMIA & operator=(SparseMIA&& otherMIA){
        this->set_dims(otherMIA.dims());

        this->setLinIdxSequence(otherMIA.linIdxSequence());
        this->setSorted(otherMIA.is_sorted());
        m_data.swap(otherMIA.m_data);
        m_indices.swap(otherMIA.m_indices);
        this->mSolveInfo=otherMIA.solveInfo();
        return *this;
    }

    //! Returns the data container
    const Data & data() const{
        return m_data;
    }
    //! Returns the index container
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

    //! Returns a raw pointer to the scalar data
    index_type* raw_index_ptr(){
        return &m_indices[0];
    }

    //! Returns a raw pointer to the scalar data
    const index_type* raw_index_ptr() const{
        return &m_indices[0];
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

    void push_back(const data_type & _data, const index_type& _index)
    {
        m_indices.push_back(_index);
        m_data.push_back(_data);
        if(this->is_sorted()&& this->size()>1){
            if(m_indices.back()<*(m_indices.end()-2))
                this->setSorted(false);
        }
    }





//    //! Converts a scalar value to data_type
//    /*!
//        \tparam from_data_type the data_type you are converting from
//    */
//    template<class from_data_type>
//    data_type convert(const from_data_type from) const{
//        using namespace boost::numeric;
//        typedef converter<data_type,from_data_type> to_mdata_type;
//        return to_mdata_type::convert(from);
//    }

    //! Returns size of raw data or the number of nonzeros
    std::size_t size() const
    {

        return m_data.size();

    }

    //! Clears data and index containers
    void clear()
    {

        m_data.clear();
        m_indices.clear();

    }

    //!resizes data and index containers. Previous references should be considered invalid
    void resize(size_t _size)
    {
        m_data.resize(_size);
        m_indices.resize(_size);
        this->mIsSorted=false;
    }

    //!reserves capacity for data and index containers. Previous references should be considered invalid
    void reserve(size_t _size)
    {

        m_data.reserve(_size);
        m_indices.reserve(_size);
    }

    //! Removes data with duplicated indices - conflicts are solved by always choosing the first data entry encountered
    void collect_duplicates()
    {
        select_first<data_type> selector;
        collect_duplicates(selector);
    }

    //! Removes data with duplicated indices - conflicts are solved by using the collector class, ie std::plus<data_type>
    template<class Collector>
    void collect_duplicates(Collector collector)
    {



        this->sort();


        auto result_idx = this->index_begin();
        auto result_data= this->data_begin();
        auto first=result_idx;
        while (++first != this->index_end())
        {
            if (*result_idx != *first){
                if(this->data_at(first)){
                    *(++result_idx)=*first;
                    *(++result_data)=this->data_at(first);
                }
            }
            else{

                *result_data=collector(*result_data,this->data_at(first));
            }
        }

        size_t diff=result_idx-this->index_begin()+1;
        resize(diff);
        this->setSorted(true);
    }


    template<typename otherDerived, typename Op,typename index_param_type>
    void merge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order);

    template<typename otherDerived, typename Op,typename index_param_type>
    void merge(DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order);


    template<typename otherDerived, typename Op>
    void merge(SparseMIABase<otherDerived> &b,const Op& op);

    template<typename otherDerived, typename Op>
    void merge(DenseMIABase<otherDerived> &b,const Op& op);

protected:

    template<class otherMIAType,typename boost::enable_if< internal::is_SparseMIA<otherMIAType>,int >::type = 0>
    void copy_other_MIA(const otherMIAType & otherMIA)
    {

        //std::cout << "copy other MIA with sparse " << std::endl;
        this->set_dims(otherMIA.dims());
        this->setLinIdxSequence(otherMIA.linIdxSequence());

        this->mSolveInfo=otherMIA.solveInfo();
        this->resize(otherMIA.size());

        auto otherIt=otherMIA.storage_begin();
        std::for_each(this->storage_begin(),this->storage_end(),[this,&otherIt](full_tuple cur_tuple){
            this->data_val(cur_tuple)=this->convert(std::get<0>(*otherIt)); //need 'this' due to bug in gcc compiler
            this->index_val(cur_tuple)=std::get<1>(*otherIt++);
        });
        this->mIsSorted=otherMIA.is_sorted();

    }

    template<class otherMIAType,typename boost::enable_if< internal::is_DenseMIA<otherMIAType>,int >::type = 0>
    void copy_other_MIA(const otherMIAType & denseMIA)
    {



        this->reset_linIdx_sequence();

        this->mSolveInfo=denseMIA.solveInfo();
        this->set_dims(denseMIA.dims());

        //count the number of nnzs
        size_t nnz=0;
        for(auto it=denseMIA.data_begin();it<denseMIA.data_end();++it)
            if(std::abs(*it)>SparseTolerance<data_type>::tolerance)
                ++nnz;

        //allocate the required data
        this->resize(nnz);


        //set m_data and m_indices

        nnz=0;
        for(auto it=denseMIA.data_begin();it<denseMIA.data_end();++it){
            if(std::abs(*it)>SparseTolerance<data_type>::tolerance){
                *(this->data_begin()+nnz)=*it;
                *(this->index_begin()+nnz++)=it-denseMIA.data_begin();
            }

        }
        this->mIsSorted=true;

    }





    //! Method to perform scanning type merge operations, eg add.
    template<typename otherDerived, typename Op,typename index_param_type>
    void scanMerge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order);

    //! Method to perform sorting type merge operations, eg add.
    template<typename otherDerived, typename Op,typename index_param_type>
    void sortMerge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order);
private:






};

template <class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type>
void  SparseMIA<T,_order>::scanMerge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order)
{

    *this=this->outside_scanMerge(b,op,index_order); //in_place merge will just require a copy anyway, so just go straight to a out-of place scan merge



//    //if b is sorted, but *this is not, then we scan based on b
//    if (b.is_sorted() && !this->is_sorted()){
//        //get the order of lhs indices in terms of rhs
//        auto lhsOrder=internal::reverseOrder(index_order);
//        //if b is also is not sorted using the default sort order, we also need to change lhsOrder to reflect this
//        lhsOrder= internal::reOrderArray(b.linIdxSequence(), lhsOrder);
//        this->sort(lhsOrder);
//    }
//    else if(this->is_sorted()&&!b.is_sorted()){
//       b.sort(internal::reOrderArray(this->mLinIdxSequence,index_order)); //change b's sort order to it matches the index order, and also *this's current sort order
//
//    } //if both *this and b are sorted, then we resort whichever mia has a smaller number of nonzeros
//    else if(this->is_sorted()&& b.is_sorted()){
//        if (b.size()<this->size()){
//            b.sort(internal::reOrderArray(this->mLinIdxSequence,index_order));
//        }
//        else{
//            //get the order of lhs indices in terms of rhs
//            auto lhsOrder=internal::reverseOrder(index_order);
//            //we also need to reorder b's sort order (which may not be {0,1,2, etc.}) using the index order
//            lhsOrder= internal::reOrderArray(b.linIdxSequence(), lhsOrder);
//            this->sort(lhsOrder);
//        }
//    }
//    else
//        throw MIAParameterException("Scan Merge should never have been called if both MIAs are unsorted");
//
//
//
//    //print_array(lhsOrder,"lhsOrder");
//    //boost::timer::cpu_timer sort_t;
//
//    //this->change_linIdx_sequence(lhsOrder);
//
//    //std::sort(this->index_begin(),this->index_end());
//    //std::sort(this->data_begin(),this->data_end());
//
//    //std::cout << "Scan sort " << boost::timer::format(sort_t.elapsed()) << std::endl;
//    //this->print();
//    size_t old_size=this->size();
//    this->resize(this->size()+b.size());
//    //boost::timer::cpu_timer merge_t;
//
//    //merge the two sorted containers to *this's containers
//    auto new_end=internal::merge_sparse_storage_containers(this->storage_begin(),this->storage_begin()+old_size,b.storage_begin(),b.storage_end(),op);
//    this->resize(new_end-this->storage_begin());

}

//both operands are assumed to be unsorted already
template <class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type>
void  SparseMIA<T,_order>::sortMerge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order)
{

    //b has less non-zeros than *this
    if(this->size()>b.size()){
        //std::cout << "this is bigger " << std::endl;
        auto b_linIdxSequence=b.linIdxSequence();

        internal::reorder_from(index_order,this->linIdxSequence(),b_linIdxSequence);

        b.change_linIdx_sequence(b_linIdxSequence); //change b's sort order to it matches the index order, and also *this's current sort order
    }
    else{
        //std::cout << "b is bigger " << std::endl;
        //get the order of lhs indices in terms of rhs
        auto lhsOrder=internal::reverseOrder(index_order);
        auto temp_linIdxSequence=this->linIdxSequence();
        internal::reorder_from(lhsOrder,b.linIdxSequence(),temp_linIdxSequence);
        this->change_linIdx_sequence(temp_linIdxSequence);
    }



    size_t old_size=this->size();
    this->resize(this->size()+b.size());
    this->mIsSorted=false;
    if(boost::is_same<std::minus<data_type>,Op>::value){
        //since stable sort is slower (inherently and also because we have to use tupleit with stable_sort), we just negate b's datatype
        //during the copy process, and then just perform addition, which doesn't need stability
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
        //don't use storage iterators, as they're likely slower than just performing two separate copy runs
        std::copy(b.data_begin(),b.data_end(),this->data_begin()+old_size);
        std::copy(b.index_begin(),b.index_end(),this->index_begin()+old_size);
        this->collect_duplicates(op);
    }




}



template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void SparseMIA<T,_order>::assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{


    //otherMIA.print();
    static_assert(_order==internal::order<otherDerived>::value,"Array specifying index shuffle must be the same size as the order of the operand array");
    auto new_dims=otherMIA.dims();
    internal::reorder_from(otherMIA.dims(),index_order,new_dims); //don't bother to re-sort, just change the linIdxSequence
    this->set_dims(new_dims);
    //print_array(index_order,"index_order");
    //print_array(otherMIA.linIdxSequence(),"otherMIA.linIdxSequence()");

    //get the order of lhs indices in terms of rhs
    auto lhsOrder=internal::reverseOrder(index_order);
    //if b is also is not sorted using the default sort order, we also need to change lhsOrder to reflect this
    internal::reorder_from(lhsOrder,otherMIA.linIdxSequence(),this->linIdxSequence());
    //print_array(this->mLinIdxSequence,"this->mLinIdxSequence");
//    if(otherMIA.linIdxSequence()!=this->linIdxSequence())
//        this->mIsSorted=false;
//    else



    this->resize(otherMIA.size());
    this->mSolveInfo=otherMIA.solveInfo();
    auto otherDataIt=otherMIA.data_begin();
    std::for_each(this->data_begin(),this->data_end(),[&](data_type & cur_data){
        cur_data=this->convert(*(otherDataIt++));
    });
    std::copy(otherMIA.index_begin(),otherMIA.index_end(),this->index_begin());
    this->mIsSorted=otherMIA.is_sorted();




}

template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void SparseMIA<T,_order>::assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{

    //std::cout << "Dense assign" << std::endl;
    static_assert(_order==internal::order<otherDerived>::value,"Array specifying index shuffle must be the same size as the order of the operand array");
    auto new_dims=this->dims();
    internal::reorder_from(otherMIA.dims(),index_order,new_dims); //shuffle the dimensions around based on index_order
    this->set_dims(new_dims);
    internal::reverseOrder(index_order,this->linIdxSequence()); //don't bother to re-sort, just change the linIdxSequence

    this->mSolveInfo=otherMIA.solveInfo();



    size_t nnz=0;
    for(auto it=otherMIA.data_begin();it<otherMIA.data_end();++it)
        if(std::abs(*it)>SparseTolerance<data_type>::tolerance)
            ++nnz;

    //allocate the required data
    this->resize(nnz);


    //set m_data and m_indices
    nnz=0;
    for(auto it=otherMIA.data_begin();it<otherMIA.data_end();++it){
        if(std::abs(*it)>SparseTolerance<data_type>::tolerance){
            *(this->data_begin()+nnz)=this->convert(*it);
            *(this->index_begin()+nnz++)=it-otherMIA.data_begin();

        }

    }
    this->setSorted(true);

}

template<class T, size_t _order>
template<typename otherDerived, typename Op>
void  SparseMIA<T,_order>::merge(DenseMIABase<otherDerived> &b,const Op& op){
    return this->merge(b,op,internal::createAscendingIndex<mOrder>());

}


template<class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type>
void SparseMIA<T,_order>::merge(DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order)
{



    this->check_merge_dims(b,index_order);

    //std::array<size_t,mOrder> converted_index_order=array_converter<size_t>::convert(index_order);
    size_t nnz=0;
    for(auto it=b.data_begin();it<b.data_end();++it)
            if(std::abs(*it)>SparseTolerance<data_type>::tolerance)
                ++nnz;
    bool negate=boost::is_same<std::minus<data_type>,Op>::value;
    auto old_size=this->size();
    this->resize(nnz+old_size);
    auto old_idx=this->index_begin()+old_size;
    auto old_data=this->data_begin()+old_size;
    for(size_t idx=0;idx<b.dimensionality();++idx){

        if(std::abs(b.atIdx(idx))>SparseTolerance<data_type>::tolerance){

            auto lhs_index=internal::sub2ind_reorder(b.ind2sub(idx),index_order,b.dims());

            if(negate)
                *old_data++=-1*b.atIdx(idx);
            else
                *old_data++=b.atIdx(idx);
            //may need to reshuffle index again, if current sort order isn't the default sort order
            *old_idx++=this->convert_from_default_sort(lhs_index);


        }
    }
    this->setSorted(false);
    if(negate)
        this->collect_duplicates(std::plus<data_type>());
    else
        this->collect_duplicates(op);



}

template<class T, size_t _order>
template<typename otherDerived, typename Op>
void  SparseMIA<T,_order>::merge(SparseMIABase<otherDerived> &b,const Op& op){
    //std::cout << "Should have got here I think " << std::endl;
    return this->merge(b,op,internal::createAscendingIndex<mOrder>());

}

template<class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type>
void  SparseMIA<T,_order>::merge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIA>::value>& index_order)
{

    //std::cout << "destructive merge " << std::endl;

    this->check_merge_dims(b,index_order);
    std::array<size_t,mOrder> converted_index_order=array_converter<size_t>::convert(index_order);



    if(b.is_sorted() || this->is_sorted()){
        //std::cout << "Performing Scan merge " <<std::endl;
        scanMerge(b,op,converted_index_order);
    }
    else{
        //std::cout << "Performing Sort merge " <<std::endl;
        sortMerge(b,op,converted_index_order);
    }
    //std::cout << "After merge " << std::endl;
    //this->print();








}


/*! @} */

}

#endif // SparseMIA_H_INCLUDED
