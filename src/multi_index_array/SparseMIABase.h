// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef SPARSEMIABASE_H_INCLUDED
#define SPARSEMIABASE_H_INCLUDED

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


#include "LibMiaException.h"
#include "Util.h"
#include "IndexUtil.h"
#include "MIA.h"



//\defgroup
namespace LibMIA
{

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */

namespace internal
{



template<class Derived>
struct data_type<SparseMIABase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<SparseMIABase<Derived> >: public index_type<Derived> {};

template<class Derived>
struct order<SparseMIABase<Derived> >: public order<Derived> {};

template<class Derived>
struct data_iterator<SparseMIABase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct const_data_iterator<SparseMIABase<Derived> >: public const_data_iterator<Derived> {};

template<class Derived>
struct index_iterator<SparseMIABase<Derived> >: public index_iterator<Derived> {};

template<class Derived>
struct const_index_iterator<SparseMIABase<Derived> >: public const_index_iterator<Derived> {};

template<class Derived>
struct storage_iterator<SparseMIABase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<SparseMIABase<Derived> >: public const_storage_iterator<Derived> {};

template<class Derived>
struct Data<SparseMIABase<Derived> >: public Data<Derived> {};

template<class Derived>
struct Indices<SparseMIABase<Derived> >: public Indices<Derived> {};

template<class Derived>
struct full_tuple<SparseMIABase<Derived> >: public full_tuple<Derived> {};

template<class Derived>
struct const_full_tuple<SparseMIABase<Derived> >: public const_full_tuple<Derived> {};

}



//!  Base class for sparse multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  sparse multi-index array types. Provides operations and functions common to all sparse
  multi-index arrays.

  \tparam Derived   should only be a sparse multi-index array class.
*/
template <class Derived>
class SparseMIABase: public MIA<SparseMIABase<Derived > >
{
public:


    typedef typename internal::data_type<SparseMIABase>::type data_type;
    typedef typename internal::index_type<SparseMIABase>::type index_type;
    typedef typename internal::Data<SparseMIABase>::type Data;
    typedef typename internal::data_iterator<SparseMIABase>::type data_iterator;
    typedef typename internal::const_data_iterator<SparseMIABase>::type const_data_iterator;
    typedef typename internal::index_iterator<SparseMIABase>::type index_iterator;
    typedef typename internal::const_index_iterator<SparseMIABase>::type const_index_iterator;
    typedef typename internal::storage_iterator<SparseMIABase>::type storage_iterator;
    typedef typename internal::const_storage_iterator<SparseMIABase>::type const_storage_iterator;
    typedef typename internal::full_tuple<SparseMIABase>::type full_tuple;
    typedef typename internal::const_full_tuple<SparseMIABase>::type const_full_tuple;
    constexpr static size_t mOrder=internal::order<SparseMIABase>::value;
    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }
    /** \returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

    template<typename... Dims>
    SparseMIABase(Dims... dims): MIA<SparseMIABase<Derived > >(dims...),mIsSorted(true) {
        init_sort_order();
    }


    SparseMIABase(): MIA<SparseMIABase<Derived > >(),mIsSorted(true) {
        init_sort_order();
    }


    SparseMIABase(const std::array<index_type,mOrder> &_dims): MIA<SparseMIABase<Derived > >(_dims),mIsSorted(true) {
        init_sort_order();
    }

    //! Returns scalar data at given linear index or first element if not found
    const data_type& atIdx(index_type idx) const{


        storage_iterator it=find_idx(idx);
        return it==derived().data_end()?*derived().data_begin():*it;
    }

    //! Returns scalar data at given linear index or first element if not found
    data_type& atIdx(index_type idx) {

        storage_iterator it=find_idx(idx);
        return it==derived().data_end()?*derived().data_begin():*it;
    }

    const_storage_iterator storage_begin() const
    {
        return boost::make_zip_iterator(boost::make_tuple(derived().data_begin(), derived().index_begin()));


    }


    const_storage_iterator storage_end() const
    {
        return boost::make_zip_iterator(boost::make_tuple(derived().data_end(), derived().index_end()));
    }

    storage_iterator storage_begin()
    {
        return iterators::makeTupleIterator(derived().data_begin(),derived().index_begin());

    }


    storage_iterator storage_end()
    {
        return iterators::makeTupleIterator(derived().data_end(),derived().index_end());
    }

    data_iterator data_begin()
    {


        return derived().data_begin();
    }

    data_iterator data_end()
    {


        return derived().data_end();
    }

    const_data_iterator data_begin() const
    {


        return derived().data_begin();
    }

    const_data_iterator data_end() const
    {


        return derived().data_end();
    }

    index_iterator index_begin()
    {


        return derived().index_begin();
    }

    index_iterator index_end()
    {


        return derived().index_end();
    }

    const_index_iterator index_begin() const
    {


        return derived().index_begin();
    }

    const_index_iterator index_end() const
    {


        return derived().index_end();
    }


    //!resizes data and index containers. Previous references should be considered invalid
    void resize(size_t _size)
    {
        derived().m_data.resize(_size);
        derived().m_indices.resize(_size);
    }

    bool is_sorted() const{
        return mIsSorted;
    }

    void sort(const std::array<size_t,mOrder> & _sort_order)
    {
        change_sort_order(_sort_order);
        sort();
    }

    void sort()
    {


        if(!mIsSorted)
        {
            std::sort(storage_begin(),storage_end(),[] (const full_tuple& left,const full_tuple& right)
            {
                return boost::get<1>(left)<boost::get<1>(right);
            } );
        }


    }

    void print() const
    {
        std::cout << "Index\t Data" << std::endl;
        for(auto it=this->storage_begin();it<this->storage_end();++it){
            std::cout << index_val(*it) << "\t " << data_val(*it) << std::endl;
        }
    }

    template<class index_param_type>
    void change_sort_order(const std::array<index_param_type,mOrder> & _sort_order)
    {
        static_assert(internal::check_index_compatibility<size_t,index_param_type>::type::value,"Must use an array convertable to index_type");
        if(_sort_order==mSortOrder)
            return;
        print_array(mSortOrder,"mSortOrder");
        for(auto& it: derived().m_indices){
            it=this->sub2ind(this->ind2sub(it,mSortOrder),_sort_order);
        }
        mSortOrder=_sort_order;
        mIsSorted=false;
    }

    void reset_sort_order(){
        change_sort_order(mDefaultSortOrder);
    }

    //!converts a linear index calculated using mSortOrder to a linear index calculated using mDefaultSortOrder
    index_type convert_to_default_sort(const index_type idx) const
    {
        return this->sub2ind(this->ind2sub(idx,mSortOrder),mDefaultSortOrder);
    }



    //!  Sets MIA index data to uniformly distributed random values.
    /*!
    May cause duplicates

    */
    void rand_indices(){
        using namespace boost::numeric;
        boost::uniform_real<> uni_dist(0,this->m_dimensionality-1);
        boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(gen, uni_dist);
        typedef converter<index_type,boost::uniform_real<>::result_type,conversion_traits<index_type,boost::uniform_real<>::result_type>,def_overflow_handler,RoundEven<boost::uniform_real<>::result_type>> to_mdata_type;
        for (auto i=derived().index_begin();i<derived().index_end();++i){
            *i=to_mdata_type::convert(uni());
        }
        mIsSorted=false;

    }

    //! Removes data with duplicated indices - conflicts are solved by always choosing the first data entry encountered
    void collect_duplicates()
    {
        select_first<data_type> selector;
        collect_duplicates(selector);
    }

    void setSorted(bool isSorted){
        mIsSorted=isSorted;
    }

    //! Removes data with duplicated indices - conflicts are solved by using the collector class, ie std::plus<data_type>
    template<class Collector>
    void collect_duplicates(Collector collector)
    {


        sort();
        auto result = storage_begin();
        auto first=result;
        while (++first != storage_end())
        {
            if (!(index_val(*result) == index_val(*first))){
                *(++result)=*first;
            }
            else{

                data_val(*result)=collector(data_val(*result),data_val(*first));
            }
        }
        ++result;
        size_t diff=result-storage_begin();
        resize(diff);
    }

    template<class otherDerived>
    bool operator==(const SparseMIABase<otherDerived>& otherMIA) const
    {

        bool passed=std::equal(this->index_begin(),this->index_end(),otherMIA.index_begin());
        if(!passed)
            return false;
        return std::equal(this->data_begin(),this->data_end(),otherMIA.data_begin());

    }

    template<class otherDerived>
    bool operator==(const DenseMIABase<otherDerived>& otherMIA) const
    {


        auto it=this->storage_begin();
        if (otherMIA.atIdx(convert_to_default_sort(index_val(*it)))!=data_val(*it))
            return false;
        for(index_type idx=0;idx<index_val(*(it));idx++)
            if (std::abs(otherMIA.atIdx(convert_to_default_sort(idx)))>SparseTolerance<data_type>::tolerance)
                return false;


        for(it=this->storage_begin()+1;it<this->storage_end();++it){
            if (otherMIA.atIdx(convert_to_default_sort(index_val(*it)))!=data_val(*it))
                return false;
            for(auto idx=index_val(*(it-1))+1;idx<index_val(*(it));idx++)
                if (std::abs(otherMIA.atIdx(convert_to_default_sort(idx)))>SparseTolerance<data_type>::tolerance)
                    return false;

        }

        for(index_type idx=*(this->index_end()-1)+1;idx<this->m_dimensionality;idx++)
            if (std::abs(otherMIA.atIdx(convert_to_default_sort(idx)))>SparseTolerance<data_type>::tolerance)
                return false;

        return true;

    }

    template<class otherDerived>
    bool operator!=(const MIA<otherDerived>& otherMIA) const
    {

        return !(*this==otherMIA.derived());

    }

    size_t size() const
    {
        return derived().size();
    }






protected:

    //!keeps track of the order of dims used to calculate linear indices
    std::array<size_t,mOrder> mSortOrder;
    //!keeps track of the order of dims used to calculate linear indices
    std::array<size_t,mOrder> mDefaultSortOrder;
    //!keeps track of whether SparseMIA is sorted or not
    bool mIsSorted;

    storage_iterator & find_idx(const index_type idx) const
    {
        return std::lower_bound(storage_begin(),storage_end(),idx,[] (const full_tuple& _tuple,const index_type idx)
        {
            return boost::get<1>(_tuple)<idx;
        } );

    }




    void init_sort_order(){
        for(size_t i=0;i<mSortOrder.size();++i)
            mSortOrder[i]=i;
        mDefaultSortOrder=mSortOrder;
    }

    void set_sort_order(const std::array<size_t,mOrder> & _sort_order){
        mSortOrder=_sort_order;
    }



    index_type& index_val(const full_tuple & a)
    {
        return boost::get<1>(a);

    }

    const index_type& index_val(const const_full_tuple & a) const
    {
        return boost::get<1>(a);

    }


    data_type& data_val(const full_tuple & a)
    {
        return boost::get<0>(a);

    }
    const data_type& data_val(const const_full_tuple & a) const
    {
        return boost::get<0>(a);

    }



    std::array<index_type,mOrder> ind2sub(index_type idx,const std::array<size_t,mOrder>& dim_order) const{
        //print_array(internal::ind2sub(idx, this->dims(),dim_order),"ind2sub");
        return internal::ind2sub(idx, this->dims(),dim_order);
    }

    index_type sub2ind(const std::array<index_type,mOrder> & indices,const std::array<size_t,mOrder>& dim_order) const{

        return internal::sub2ind_reorder(indices, dim_order,this->dims());
    }

    template<typename otherDerived, typename Op,typename index_param_type>
    void merge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order);

private:

    friend class MIA<SparseMIABase>;

};

template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type>
void  SparseMIABase<Derived>::merge(const SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order)
{

    this->check_merge_dims(b,index_order);
    std::array<size_t,mOrder> converted_index_order=array_converter<size_t>::convert(index_order);



    if(b.is_sorted()){
        //std::cout << "Performing Scan merge " <<std::endl;
        derived().scanMerge(b,op,converted_index_order);
    }
    else{
        derived().sortMerge(b,op,converted_index_order);
    }







}



/*! @} */

}

#endif // SPARSEMIABASE_H_INCLUDED
