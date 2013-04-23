// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef SPARSELATTICEBASE_H
#define SPARSELATTICEBASE_H
#include <algorithm>
#include <vector>
#include <type_traits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
//#include <functional>

//#include <boost/mpl/assert.hpp>
//#include <boost/mpl/if.hpp>
//#include <boost/mpl/bool.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/numeric/conversion/converter.hpp>



#include "MIAConfig.h"
#include <Eigen/Sparse>
//if we're using sparse_solve then we must include SuperLU support
#ifdef LIBMIA_USE_SPARSE_SOLVE
#include <Eigen/SuperLUSupport>
#endif


#include "Lattice.h"
#include "LibMIAAlgorithm.h"
#include "Util.h"
#include "tupleit.hh"
//#define LM_CSC_TIMES 1 //define to perform old compressed column lattice multiplication - should be off
//#define LM_COLUMN_SEARCH 1 //define to perform column search during mult_scatter operation (not as efficient)
namespace LibMIA
{

namespace internal
{


template<class Derived>
struct data_type<SparseLatticeBase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<SparseLatticeBase<Derived> >: public index_type<Derived> {};

template<class Derived>
struct index_iterator<SparseLatticeBase<Derived> >: public index_iterator<Derived> {};

template<class Derived>
struct data_iterator<SparseLatticeBase<Derived> >: public data_iterator<Derived> {};

}

/** \addtogroup lattice Lattice Classes
 *  @{
*/

//!  Lattice class for sparse data.
/*!
  Supports addition, multiplication, and solution of, possible over-determined, systems of
  linear equations.


  \tparam T   the datatype of individual elements.
  \tparam managed if false, lattice is mapped from existing data. Upon destruction, lattice will not
                    free data. Default is true.
*/
template <class Derived>
class SparseLatticeBase: public Lattice<SparseLatticeBase<Derived> >
{
public:


    typedef typename internal::data_type<Derived>::type data_type;
    typedef typename LibMIA::internal::index_type<Derived>::type index_type;
    typedef typename LibMIA::internal::Data<Derived>::type Data;
    typedef typename LibMIA::internal::Indices<Derived>::type Indices;
    typedef typename LibMIA::internal::full_iterator_tuple<Derived>::type full_iterator_tuple;
    typedef typename LibMIA::internal::const_full_iterator_tuple<Derived>::type const_full_iterator_tuple;
    typedef typename LibMIA::internal::full_tuple<Derived>::type full_tuple;
    typedef typename LibMIA::internal::const_full_tuple<Derived>::type const_full_tuple;
    typedef typename LibMIA::internal::storage_iterator<Derived>::type storage_iterator;
    typedef typename LibMIA::internal::const_storage_iterator<Derived>::type const_storage_iterator;
    typedef typename LibMIA::internal::index_iterator<Derived>::type index_iterator;
    typedef typename LibMIA::internal::const_index_iterator<Derived>::type const_index_iterator;
    typedef typename LibMIA::internal::data_iterator<Derived>::type data_iterator;
    typedef typename LibMIA::internal::const_data_iterator<Derived>::type const_data_iterator;
    typedef typename Eigen::SparseMatrix<data_type> SparseMatrix_cm;
    typedef typename Eigen::SparseMatrix<data_type,Eigen::RowMajor> SparseMatrix_rm;
    typedef typename Eigen::MappedSparseMatrix<data_type,Eigen::ColMajor,index_type> MappedSparseMatrix_cm;
    typedef typename Eigen::MappedSparseMatrix<data_type,Eigen::RowMajor,index_type> MappedSparseMatrix_rm;

    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }
    /** \returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

    //only enable for operands that have the same index_type
    template <class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type operator*(SparseLatticeBase<otherDerived> &b);


    //only enable for operands that have the same index_type
    template <class otherDerived>
    typename SparseSolveReturnType<Derived,otherDerived>::type solve(SparseLatticeBase<otherDerived> &b);

    data_type operator()(index_type _row, index_type _column, index_type _tab) const;





//
//

//


    void sparse_init(bool _is_sorted, bool _sort_order)
    {

        m_is_sorted=_is_sorted;
        m_sort_order=_sort_order;


    }

    template<class otherDerived>
    typename SparseMergeReturnType<Derived,otherDerived>::type operator+(SparseLatticeBase<otherDerived> &b)
    {
        typedef typename SparseMergeReturnType<Derived,otherDerived>::type RType;
        RType c(*this);
        c+=b;
        return c;

    }

    template<class otherDerived>
    typename SparseMergeReturnType<Derived,otherDerived>::type operator-(SparseLatticeBase<otherDerived> &b)
    {
        typedef typename SparseMergeReturnType<Derived,otherDerived>::type RType;
        RType c(*this);
        c-=b;
        return c;

    }


    void print()
    {





        sort(m_sort_order);
        std::cout << "Val\t" ;
        std::cout << "Row\t"<< "Column\t" << "Tab\n";
        for (storage_iterator i=storage_begin(); i<storage_end(); i++)
        {

            full_tuple temp=*i;
            //T temp2=std::get<0>(temp);

            std::cout << std::get<0>(temp) <<"\t" ;
            std::cout << row(index_val(temp)) << "\t"<< column(index_val(temp)) << "\t" << tab(index_val(temp)) <<"\n";
        }
    }

    std::size_t size() const
    {
        return derived().size();

    }
    bool idx_less(index_type a, index_type b) const
    {


        return a<b;

    }

    bool idx_less_rowmajor(index_type a, index_type b) const
    {


        return (column(a)+row(a)*this->width()+tab(a)*this->width()*this->height())<(column(b)+row(b)*this->width()+tab(b)*this->width()*this->height());

    }


    void sort(bool sort_order=ColumnMajor)
    {





        storage_iterator temp=storage_begin();

        if (m_sort_order!=sort_order || !is_sorted())
        {

            if (sort_order==ColumnMajor)
                internal::Introsort(this->index_begin(),this->index_end(),std::less<index_type>(),internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));

            else{
                std::function<bool(const index_type&,const index_type&)> comparator=[this](const index_type & lhs, const index_type & rhs){
                    return this->idx_less_rowmajor(lhs,rhs);
                };
                internal::Introsort(this->index_begin(),this->index_end(),comparator,internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));

            }

            m_is_sorted=true;
            m_sort_order=sort_order;
        }


    }


    //!  Sets Lattice index data to uniformly distributed random values.
    /*!
    May cause duplicates

    */
    void rand_indices(){
        using namespace boost::numeric;

        boost::uniform_real<> uni_dist(0,this->dimensionality()-1);
        boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(gen, uni_dist);
        typedef converter<index_type,boost::uniform_real<>::result_type,conversion_traits<index_type,boost::uniform_real<>::result_type>,def_overflow_handler,RoundEven<boost::uniform_real<>::result_type>> to_mdata_type;
        for (auto i=derived().index_begin();i<derived().index_end();++i){
            *i=to_mdata_type::convert(uni());
        }
        m_is_sorted=false;

    }




    index_type row(index_type lin_index) const
    {

        return lin_index%this->height();

    }

    index_type column(index_type lin_index) const
    {
        lin_index/=this->height();
        return lin_index%this->width();

    }

    index_type tab(index_type lin_index) const
    {
        lin_index/=this->height();
        lin_index/=this->width();
        return lin_index%this->depth();

    }




    bool is_sorted() const
    {

        return m_is_sorted;
    }
    bool sort_order() const
    {

        return m_sort_order;
    }

    bool tab_compare(index_type lin_index1, index_type lin_index2) const
    {
        return tab(lin_index1)<tab(lin_index2);

    }
    bool column_compare(index_type lin_index1, index_type lin_index2) const
    {
        return column(lin_index1)<column(lin_index2);

    }
    bool row_compare(index_type lin_index1, index_type lin_index2) const
    {
        return row(lin_index1)<row(lin_index2);

    }


    index_type& index_val(full_tuple a)
    {
        return std::get<1>(a);

    }

    const index_type& index_val(const_full_tuple a) const
    {
        return std::get<1>(a);

    }



    data_type& data_val(full_tuple a)
    {
        return std::get<0>(a);

    }
    const data_type& data_val(const_full_tuple a) const
    {
        return std::get<0>(a);

    }

    const data_type& data_at(size_t nnz_index) const
    {
        return *(data_begin()+nnz_index);

    }

    data_type& data_at(size_t nnz_index)
    {
        return *(data_begin()+nnz_index);

    }

    const index_type& index_at(size_t nnz_index) const
    {
        return *(index_begin()+nnz_index);

    }

    index_type& index_at(size_t nnz_index)
    {
        return *(index_begin()+nnz_index);

    }




    index_iterator index_begin()
    {
        return derived().index_begin();

    }

    const_index_iterator index_begin() const
    {
        return derived().index_begin();

    }

    index_iterator index_end()
    {
        return derived().index_end();

    }

    const_index_iterator index_end() const
    {
        return derived().index_end();

    }

    data_iterator data_begin()
    {
        return derived().data_begin();

    }

    const_data_iterator data_begin() const
    {
        return derived().data_begin();

    }

    data_iterator data_end()
    {
        return derived().data_end();

    }

    const_data_iterator data_end() const
    {
        return derived().data_end();

    }



    storage_iterator storage_begin();

    const_storage_iterator storage_begin() const;

    storage_iterator storage_end();

    const_storage_iterator storage_end() const;

    //non constant b/c a sort may be involved
    template<class otherDerived>
    bool operator==(const DenseLatticeBase<otherDerived>& otherLat);

    //non constant b/c a sort may be involved
    template<class otherDerived>
    bool fuzzy_equals(const DenseLatticeBase<otherDerived>& otherLat,double precision);

    //non constant b/c a sort may be involved
    template<class otherDerived>
    bool operator!=(const DenseLatticeBase<otherDerived>& otherLat){
        return !(*this==otherLat);
    }

    //!  Utility function for mapping tabs to column-major matrices to use in lattice mulitplication
    /*!
    \param[in] row_map A sorted vector mapping compressed indices to full row indices of the current tab
    \param[in] inner_indices A sorted vector mapping shared column indices to full column indices of the current tab
    \param[in] t_begin t_end Storage iterator to the nonzeros of the current tab
    \param[out] mat Resulting CSC matrix. The CSC matrix should already have allocated vectors for data and indices.
    */

    void to_matrix(const std::vector<index_type> &row_map,const std::vector<index_type> &inner_indices,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_cm& mat) const;


    //!  Utility function for mapping tabs to row-major matrices to use in lattice mulitplication. See to_matrix for more information.

    void to_matrix_rowmajor(const std::vector<index_type> &inner_indices,const std::vector<index_type> &column_map,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_rm& mat) const;

    index_iterator find_tab_start_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag=true);
    index_iterator find_tab_end_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag=true);

protected:

    template<class b_index_type, class b_data_type,class ret_index_type,class super_data_type>
    index_iterator mult_scatter(std::vector<size_t>::iterator & a_column_idx_begin,std::vector<size_t>::iterator a_column_idx_end,index_iterator a_cur_it,index_iterator a_end, b_index_type b_row,b_index_type b_column,b_data_type beta,
                                                        std::vector<size_t>& row_marker,std::vector<super_data_type> & data_collector,std::vector<ret_index_type> & c_indices);


    template<class otherDerived, class BinaryPredicate>
    bool compare_with_dense(const DenseLatticeBase<otherDerived>& otherLat,BinaryPredicate predicate);

    index_type full2lin_index(index_type _row, index_type _column, index_type _tab) const;
    bool m_is_sorted;
    bool m_sort_order;


};

template<typename Derived>
inline typename SparseLatticeBase<Derived>::index_type
SparseLatticeBase<Derived>::full2lin_index(index_type _row, index_type _column, index_type _tab) const
{

    return _row+_column*this->height()+_tab*this->height()*this->width();

}

template<typename Derived>
inline typename SparseLatticeBase<Derived>::data_type
SparseLatticeBase<Derived>::operator()(index_type _row, index_type _column, index_type _tab) const
{
    sort(ColumnMajor);
    const_index_iterator found=std::find(index_begin(),index_end(),full2lin_index( _row,  _column,  _tab));

    if (found==index_end())
        return 0;
    else
        return data_at(found-index_begin());


}

template<typename Derived>
inline typename SparseLatticeBase<Derived>::storage_iterator
SparseLatticeBase<Derived>::storage_begin()
{


    return iterators::makeTupleIterator(data_begin(),index_begin());



}


template<typename Derived>
inline typename SparseLatticeBase<Derived>::const_storage_iterator
SparseLatticeBase<Derived>::storage_begin() const
{


    return iterators::makeTupleIterator(data_begin(),index_begin());



}


template<typename Derived>
inline typename SparseLatticeBase<Derived>::storage_iterator
SparseLatticeBase<Derived>::storage_end()
{


    return iterators::makeTupleIterator(data_end(),index_end());



}


template<typename Derived>
inline typename SparseLatticeBase<Derived>::const_storage_iterator
SparseLatticeBase<Derived>::storage_end() const
{



    return iterators::makeTupleIterator(data_end(),index_end());


}


template<class Derived>
template<class otherDerived>
bool SparseLatticeBase<Derived>::operator==(const DenseLatticeBase<otherDerived>& otherLat)
{
    typedef typename DenseLatticeBase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[](data_type a,other_data_type b){
        return a==b;
    };
    return compare_with_dense(otherLat,pred);
}

template<class Derived>
template<class otherDerived>
bool SparseLatticeBase<Derived>::fuzzy_equals(const DenseLatticeBase<otherDerived>& otherLat,double precision)
{
    typedef typename DenseLatticeBase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[precision](data_type a,other_data_type b){
        return std::abs(a-b)<precision;
    };
    return compare_with_dense(otherLat,pred);
}

template<class Derived>
template<class otherDerived, class BinaryPredicate>
bool SparseLatticeBase<Derived>::compare_with_dense(const DenseLatticeBase<otherDerived>& otherLat,BinaryPredicate predicate)
{
    if(this->height()!=otherLat.height() || this->width()!=otherLat.width() || this->depth()!=otherLat.depth())
        return false;


    if (!this->size())
    {
        for(auto it=otherLat.data_begin(); it<otherLat.data_end(); ++it)
            if(std::abs(*it)>SparseTolerance<data_type>::tolerance)
                return false;
        return true;
    }
    else{
        this->sort(ColumnMajor);
        auto it=this->storage_begin();
        if (!predicate(otherLat.atIdx(index_val(*it)),data_val(*it)))
            return false;

        for(index_type idx=0; idx<index_val(*(it)); idx++)
            if (std::abs(otherLat.atIdx(idx))>SparseTolerance<data_type>::tolerance)
                return false;


        for(it=this->storage_begin()+1; it<this->storage_end(); ++it)
        {
            if (!predicate(otherLat.atIdx(index_val(*it)),data_val(*it)))
            {
                //std::cout << "Trigered " << index_val(*it) << " " << convert_to_default_sort(index_val(*it)) << " " << data_val(*it) << " " << otherMIA.atIdx(convert_to_default_sort(index_val(*it))) << std::endl;

                return false;
            }
            for(auto idx=index_val(*(it-1))+1; idx<index_val(*(it)); idx++)
                if (std::abs(otherLat.atIdx(idx))>SparseTolerance<data_type>::tolerance)
                    return false;

        }

        for(index_type idx=*(this->index_end()-1)+1; idx<this->height()*this->width()*this->depth(); idx++)
            if (std::abs(otherLat.atIdx(idx))>SparseTolerance<data_type>::tolerance)
                return false;

        return true;
    }

}

template <class Derived>
auto SparseLatticeBase<Derived>::find_tab_start_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag)->index_iterator
{
    if(search_flag)
        start_it= std::lower_bound(start_it,end_it,_tab,[this](index_type lhs, index_type rhs){
                return this->tab(lhs)<rhs;
            });
    else{
        while(start_it<end_it && this->tab(*start_it)<_tab)
            start_it++;
    }
    return start_it;

}

template <class Derived>
auto SparseLatticeBase<Derived>::find_tab_end_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag)->index_iterator
{
    if(search_flag)
        start_it= std::upper_bound(start_it,end_it,_tab,[this](index_type lhs, index_type rhs){
                return lhs<this->tab(rhs);
            });
    else{
        while(start_it<end_it && this->tab(*start_it)<=_tab)
            start_it++;
    }
    return start_it;

}
#ifndef LM_CSC_TIMES
template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::operator*(SparseLatticeBase<otherDerived> &b){

    this->check_mult_dims(b);
    typedef typename ScalarPromoteType<Derived,otherDerived>::type super_data_type;

    typedef typename SparseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions
    //iterators to the current tab start and end indexes
    auto a_temp_begin=this->index_begin(), a_temp_end=this->index_begin();
    auto a_index_end=this->index_end();
    auto b_temp_begin=b.index_begin(),b_temp_end=b.index_begin();
    auto b_index_end=b.index_end();
    //maps the rows of each tab of A from compressed indices to full indices
    std::vector<index_type> a_row_map;
    //stores the starting location of each unique column in each tab of A
    std::vector<size_t> a_column_idx;
    //workspace vectors (similar to what's done when CSC matrices are multipled - see Direct Methods for Sparse Linear Systems
    std::vector<size_t> row_marker;
    std::vector<super_data_type> data_collector;
    index_type old_m=this->height();
    index_type cur_tab;
    //initial estimate of size of C
    typename RType::Indices c_indices;
    c_indices.reserve(this->size()+b.size());
    typename RType::Data c_data;
    c_data.reserve(this->size()+b.size());


    this->sort();
    b.sort();

    //determine whether we want to binary search or just scan through elements
    bool a_search_flag=false, b_search_flag=false;
    if(this->depth()*std::log(this->size())<this->size())
        a_search_flag=true;
    if(this->depth()*std::log(b.size())<b.size())
        b_search_flag=true;

    while(a_temp_begin<a_index_end && b_temp_begin<b_index_end){
        //if we're at the same tab, no work
        if(this->tab(*a_temp_begin)==b.tab(*b_temp_begin)){
            cur_tab=this->tab(*a_temp_begin);
        }
        else if (this->tab(*a_temp_begin)<b.tab(*b_temp_begin)){ //if a's tab is less than b's tab - we need to try to find b's tab in a
            cur_tab=b.tab(*b_temp_begin);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
            std::cout << "A-- Tab " << this->tab(*a_temp_begin) << " is less than b tab: " << b.tab(*b_temp_begin) <<std::endl;
#endif
            a_temp_begin=this->find_tab_start_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
            std::cout   << " so searched for it and got index " << a_temp_begin-this->index_begin() << std::endl;
#endif
            if(a_temp_begin==a_index_end) //no tab in A is greater than or equal to B's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in A equal to B's current tab, but we found one greater than it - so now we need to try to find a matching tab in b
            if(this->tab(*a_temp_begin)!=b.tab(*b_temp_begin))
                continue;
        }
        else{
            cur_tab=this->tab(*a_temp_begin);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
            std::cout << "B-- Tab " << b.tab(*b_temp_begin) << " is less than A tab: " << this->tab(*a_temp_begin) <<std::endl;
#endif
            b_temp_begin=b.find_tab_start_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
            std::cout   << " so searched for it and got index " << b_temp_begin-b.index_begin() << std::endl;
#endif
            if(b_temp_begin==b_index_end) //no tab in B is greater than or equal to A's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in B equal to A's current tab, but we found one greater than it - so now we need to try to find a matching tab in B
            if(this->tab(*a_temp_begin)!=b.tab(*b_temp_begin))
                continue;
        }

        //find the end of the current tab
        a_temp_end=this->find_tab_end_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
        std::cout << "A-- Tab " << cur_tab << " begin: " << a_temp_begin-this->index_begin() << " end: " << a_temp_end-this->index_begin() << std::endl;
#endif
        b_temp_end=b.find_tab_end_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
        std::cout << "B-- Tab " << cur_tab << " begin: " << b_temp_begin-b.index_begin() << " end: " << b_temp_end-b.index_begin() << std::endl;
#endif

        std::sort(a_temp_begin,a_temp_end,[this](index_type lhs,index_type rhs)
        {
            return this->row(lhs)<this->row(rhs);
        });
        //find the number of unique rows
        size_t unique_counter=1;
        for(auto it=a_temp_begin+1; it<a_temp_end; ++it)
        {
            if (this->row(*it)!=this->row(*(it-1)))
                unique_counter++;
        }
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
        std::cout << cur_tab << " unique_counter " << unique_counter << std::endl;
#endif
        size_t new_m=unique_counter;
        //resize our row map
        a_row_map.resize(new_m);
        a_row_map.assign(new_m,0);
        //also resize our workspace
        row_marker.resize(new_m);
        row_marker.assign(new_m,0);
        data_collector.resize(new_m);
        unique_counter=0;
        //create a row map from compressed rows to full rows, and also map the current tab's row indices to the compressed form
        index_type old_row=this->row(*a_temp_begin);
        a_row_map[0]=old_row;
        *a_temp_begin=new_m*(this->column(*a_temp_begin)+this->tab(*a_temp_begin)*this->width());
        for(auto it=a_temp_begin+1; it<a_temp_end; ++it)
        {
            if (old_row!=this->row(*it))
            {
                old_row=this->row(*it);
                unique_counter++;
                a_row_map[unique_counter]=old_row;
            }
            *it=unique_counter+new_m*(this->column(*it)+this->tab(*it)*this->width());

        }
        //temporarily set A's height to compressed number of rows
        this->m_height=new_m;
        //sort back to column major
        std::sort(a_temp_begin,a_temp_end);
#ifdef LM_COLUMN_SEARCH
        //find the number of unique columns
        unique_counter=1;
        for(auto it=a_temp_begin+1; it<a_temp_end; ++it)
        {
            if (this->column(*it)!=this->column(*(it-1)))
                unique_counter++;
        }

        //resize our column idx marker
        a_column_idx.resize(unique_counter);
        unique_counter=1;
        //store the location of the start of the first column
        a_column_idx[0]=a_temp_begin-this->index_begin();
        for(auto it=a_temp_begin+1; it<a_temp_end; ++it)
        {
            if (this->column(*it)!=this->column(*(it-1)))
            {
                //store the where in the index array the new column starts
                a_column_idx[unique_counter++]=it-this->index_begin();
            }
        }
#endif
        //iterate through every element of b
        b_index_type cur_column;
        auto cur_b=b_temp_begin;
        while(cur_b<b_temp_end)
        {
            cur_column=b.column(*cur_b);
            //for each column of b, we start at the beginning of A
            index_iterator a_cur_it=a_temp_begin;
            auto a_cur_column_idx=a_column_idx.begin();
            size_t old_c_size=c_indices.size();
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
            std::cout << "Tab " << cur_tab << " Column " << cur_column << ": old_c_size " << old_c_size << std::endl;
#endif
            while(cur_b<b_temp_end && b.column(*cur_b)==cur_column)
            {
                //should be the

                a_cur_it=mult_scatter(a_cur_column_idx,a_column_idx.end(),a_cur_it,a_temp_end,b.row(*cur_b),cur_column,b.data_at(cur_b-b.index_begin()),
                                          row_marker,data_collector,c_indices);
                cur_b++;
                //if we've reached the end of A's current tab, we're done looking at b's current column
                if(a_cur_it==a_temp_end)
                    break;
            }
            //if we've added to c's indices in the current column of B, then clean it up and add to c's data
            //note mult_scatter just pushes the rows of c to c_indices
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
            std::cout << "Tab " << cur_tab << " Column " << cur_column << ": new_c_size " << c_indices.size() << std::endl;
#endif
            if(c_indices.size()-old_c_size)
            {
                std::sort(c_indices.begin()+old_c_size,c_indices.end()); //scatter doesn't put the rows in sorted order
                c_data.resize(c_indices.size());
                for(auto it=c_indices.begin()+old_c_size; it<c_indices.end(); ++it)
                {
                    c_data[it-c_indices.begin()]=data_collector[(size_t)*it];
                    *it=a_row_map[(size_t)*it]+old_m*(cur_column+cur_tab*b.width()); //decompresses the index
                }
            }

        }
        //decompress A's indices
        for(auto it=a_temp_begin; it<a_temp_end; ++it)
        {
            *it=a_row_map[this->row(*it)]+old_m*(this->column(*it)+this->width()*cur_tab);
        }
        //reset A's height to the original height
        this->m_height=old_m;





        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;

    }

    return RType(std::move(c_data),std::move(c_indices),this->height(),b.width(),this->depth());

}
#endif
template<class Derived>
template<class b_index_type, class b_data_type,class ret_index_type,class super_data_type>
auto SparseLatticeBase<Derived>::mult_scatter(std::vector<size_t>::iterator & a_column_idx_begin,std::vector<size_t>::iterator a_column_idx_end,index_iterator a_cur_it,index_iterator a_end, b_index_type b_row,b_index_type b_column,b_data_type beta,
                                                        std::vector<size_t>& row_marker,std::vector<super_data_type> & data_collector,std::vector<ret_index_type> & c_indices)->index_iterator
{
    if (this->column(*a_cur_it)!=b_row){
        #ifdef LM_COLUMN_SEARCH
        //search the column map for the start of the column equalling b's row
        a_column_idx_begin=std::lower_bound(a_column_idx_begin,a_column_idx_end,b_row,[this](size_t lhs,b_index_type rhs){
            return this->column(this->index_at(lhs))<rhs;
        });
        //if the column we're looking for is greater than all remaining columns, just return the end of the current tab
        if (a_column_idx_begin==a_column_idx_end)
            a_cur_it= a_end;
        else{
            //otherwise, we've found a column equal to or greater than we're looking for, update our current iterator
            a_cur_it=this->index_begin()+(*a_column_idx_begin);
        }



        #else

        //search the column map for the start of the column equalling b's row
        //std::cout << "Search Range " << a_end-a_cur_it << std::endl;
        a_cur_it=std::lower_bound(a_cur_it,a_end,b_row,[this](index_type lhs,b_index_type rhs){
            return this->column(lhs)<rhs;
        });
        //std::cout << "Searching for column equal to " << b_row << " and returned " << a_cur_it-this->index_begin() << " to get a column of " << this->column(*a_cur_it) << std::endl;
        #endif



        //if we didn't find the column, but found one greater than it or reached the end, update the current index_iterator of A
        if(a_cur_it==a_end || this->column(*a_cur_it)!=b_row)
            return a_cur_it;

    }
    //std::cout << "Index of column at " << b_row << " is " << a_cur_it-this->index_begin() << std::endl;
    size_t cur_a_row;
    while(a_cur_it<a_end && this->column(*a_cur_it)==b_row){
        cur_a_row=(size_t)this->row(*a_cur_it); //here we've assumed a's rows have already been truncated
        if (row_marker[cur_a_row]<b_column+1){
            row_marker[cur_a_row]=b_column+1;
            c_indices.push_back(cur_a_row); //push back the compressed row (this will have to be remapped afterwards)
            //std::cout << "Beta: " << beta << " a_data " << this->data_at(a_cur_it-this->index_begin()) << " index " << a_cur_it-this->index_begin() << " result " << beta*this->data_at(a_cur_it-this->index_begin()) << std::endl;
            data_collector[cur_a_row]=beta*this->data_at(a_cur_it-this->index_begin());
        }
        else{
            //std::cout << "Beta: " << beta << " a_data " << this->data_at(a_cur_it-this->index_begin()) << " index " << a_cur_it-this->index_begin() << " result " << beta*this->data_at(a_cur_it-this->index_begin()) << std::endl;
            data_collector[cur_a_row]+=beta*this->data_at(a_cur_it-this->index_begin());
        }

        a_cur_it++;

    }
    return a_cur_it;

}


#ifdef LM_CSC_TIMES
template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::operator*(SparseLatticeBase<otherDerived> &b)
{

//    std::cout << "A Lattice in mult " <<std::endl;
//    this->print();
//    std::cout << "B Lattice in mult " <<std::endl;
//    b.print();
    this->check_mult_dims(b);


    typedef typename ScalarPromoteType<Derived,otherDerived>::type super_data_type;

    typedef typename SparseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions

    otherDerived b_d=b.derived();
    Derived a_d=this->derived();
    sort(ColumnMajor); //tab/column major for A
    bool old_sort_order=b.sort_order();
    b.sort(RowMajor); //tab/row major for B

    //iterators for indices and data
    storage_iterator a_begin=storage_begin();
    index_iterator a_index_begin=a_d.index_begin();
    data_iterator a_data_begin=data_begin();

    auto b_begin=b.storage_begin();
    auto b_index_begin=b.index_begin();
    auto b_data_begin=b.data_begin();

    typename RType::Indices c_indices;
    typename RType::Data c_data;

    index_iterator a_index_end=a_d.index_end();
    auto b_index_end=b.index_end();

    index_iterator a_temp_begin=a_index_begin;
    auto b_temp_begin=b_index_begin;
    index_iterator a_temp_end;
    decltype(b_temp_begin) b_temp_end;

    //iterator for row and column indices
    std::vector<index_type> a_columns;
    std::vector<index_type> a_rows;
    std::vector<b_index_type> b_columns;
    std::vector<b_index_type> b_rows;
    std::vector<index_type> inner_indices;

    typename std::vector<index_type>::iterator a_column_end;
    typename std::vector<index_type>::iterator a_row_end;
    typename std::vector<b_index_type>::iterator b_column_end;
    typename std::vector<b_index_type>::iterator b_row_end;
    typename std::vector<index_type>::iterator inner_end;

    index_type (SparseLatticeBase::*column)(index_type lin_index) const = &SparseLatticeBase::column;
    index_type (SparseLatticeBase::*row)(index_type lin_index) const = &SparseLatticeBase::row;
    b_index_type (otherDerived::*other_row)(b_index_type lin_index) const = &otherDerived::row;
    b_index_type (otherDerived::*other_column)(b_index_type lin_index) const = &otherDerived::column;
    index_type (SparseLatticeBase::*tab)(index_type lin_index) const = &SparseLatticeBase::tab;
    b_index_type (otherDerived::*tab_other)(b_index_type lin_index) const = &otherDerived::tab;
    for (int k=0; k<this->depth(); k++) //loop through every tab
    {


        //find the last occurence of current tab in both lattices
        a_temp_end=std::upper_bound(a_temp_begin,a_index_end,k,boost::bind(std::less<index_type>(), boost::lambda::_1, boost::bind(tab,this,boost::lambda::_2)));
        b_temp_end=std::upper_bound(b_temp_begin,b_index_end,k,boost::bind(std::less<b_index_type>(), boost::lambda::_1, boost::bind(tab_other,&b_d,boost::lambda::_2)));

        //if both tabs have nonzeros, then perform matrix multiplication
        if (a_temp_end!=a_temp_begin && b_temp_end!=b_temp_begin)
        {

            //get unique columns of A //***todo - just use a binary predicate in std::unique_copy, avoiding having to store needless values
            a_columns.resize(a_temp_end-a_temp_begin);
            std::transform(a_temp_begin,a_temp_end,a_columns.begin(),boost::bind(column, this,boost::lambda::_1));
            a_column_end=std::unique(a_columns.begin(),a_columns.end());
            a_columns.resize( a_column_end - a_columns.begin() );


            //get unique and sorted rows of A. This creates a map from compressed indices to actual ones //***todo - just a series of set unions to create a unique row index
            a_rows.resize(a_temp_end-a_temp_begin);
            std::transform(a_temp_begin,a_temp_end,a_rows.begin(),boost::bind(row, this,boost::lambda::_1));
            std::sort(a_rows.begin(),a_rows.end());
            a_row_end=std::unique(a_rows.begin(),a_rows.end());
            a_rows.resize( a_row_end - a_rows.begin() );



            //get unique rows of B and a map from compressed indices to actual ones
            b_rows.resize(b_temp_end-b_temp_begin);
            std::transform(b_temp_begin,b_temp_end,b_rows.begin(),boost::bind(other_row, &b_d,boost::lambda::_1));
            b_row_end=std::unique(b_rows.begin(),b_rows.end());
            b_rows.resize(b_row_end - b_rows.begin());

            //get unique and sorted columns of B. This creates a map from compressed indices to actual ones
            b_columns.resize(b_temp_end-b_temp_begin);
            std::transform(b_temp_begin,b_temp_end,b_columns.begin(),boost::bind(other_column, &b_d,boost::lambda::_1));
            std::sort(b_columns.begin(),b_columns.end());
            b_column_end=std::unique(b_columns.begin(),b_columns.end());
            b_columns.resize( b_column_end - b_columns.begin() );


            //find union of A columns and B rows
            inner_indices.resize(a_columns.size()+b_rows.size());
            inner_end=std::set_union(a_columns.begin(),a_columns.end(),b_rows.begin(),b_rows.end(),inner_indices.begin());
            inner_indices.resize(inner_end-inner_indices.begin());



            //create temporary sparse matrices for current tab and calculate matrix product
            a_columns.resize(inner_indices.size()+1);
            std::vector<index_type> rows_place_holder;
            rows_place_holder.resize(a_temp_end-a_temp_begin);
            //std::cout << "Inner size " << inner_indices.size() << " Rows size " << a_rows.size() << " data size " << a_temp_end-a_temp_begin << std::endl;
            MappedSparseMatrix_cm A(a_rows.size(),inner_indices.size(),a_temp_end-a_temp_begin,&a_columns[0],&rows_place_holder[0],&(*(a_data_begin+(a_temp_begin-a_index_begin))));

            to_matrix(a_rows,inner_indices,a_begin+(a_temp_begin-a_index_begin),a_begin+(a_temp_end-a_index_begin),A);
            //std::cout << "A Matrix in mult " <<std::endl;
            //std::cout << A <<std::endl;
            //std::cout << "Finished to _matrix" << std::endl;
            b_rows.resize(inner_indices.size()+1);
            std::vector<b_index_type> columns_place_holder;
            columns_place_holder.resize(b_temp_end-b_temp_begin);
            MappedSparseMatrix_rm B(inner_indices.size(),b_columns.size(),b_temp_end-b_temp_begin,&b_rows[0],&columns_place_holder[0],&(*(b_data_begin+(b_temp_begin-b_index_begin))));
            b.to_matrix_rowmajor(inner_indices,b_columns,b_begin+(b_temp_begin-b_index_begin),b_begin+(b_temp_end-b_index_begin),B);
            //std::cout << "Finished to _matrix rowmajor" << std::endl;

            //std::cout << "B Matrix in mult " <<std::endl;
            //std::cout << B <<std::endl;

            SparseMatrix_cm C=A*B;

            //std::cout << "C Matrix in mult " <<std::endl;
            //std::cout << C << std::endl;

            //std::cout << "About to make C tab " <<std::endl;
            //store result in two vectors - prune while doing so
            for (int j=0; j<C.outerSize(); ++j)
            {
                for (typename SparseMatrix_cm::InnerIterator it(C,j); it; ++it)
                {

                    if (std::abs(it.value())>SparseTolerance<data_type>::tolerance)
                    {
                        c_data.push_back(it.value());
                        c_indices.push_back(a_rows[it.row()]+b_columns[it.col()]*this->height()+k*this->height()*b.width()); //decompresses row and column indices
                    }

                }
            }
            //std::cout << "Finished to make C tab " <<std::endl;


        }
        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;

    }
    //std::cout << "About to sort " <<std::endl;
    //resort B in the proper precedence
    b_d.sort(old_sort_order);
    //std::cout << "Finished Mult " <<std::endl;
    //return sparse lattice using the data and indices vectors
    return RType(std::move(c_data),std::move(c_indices),this->height(),b.width(),this->depth());
}
#endif
template <class Derived>
inline void SparseLatticeBase<Derived>::to_matrix(const std::vector<index_type> &row_map, const std::vector<index_type> &inner_indices,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_cm& mat) const
{





    index_type outer_val=0;
    auto mat_rows=mat.innerIndexPtr();
    auto mat_columns=mat.outerIndexPtr();
    typename std::vector<index_type>::const_iterator cur_row;
    for(auto cur_it=inner_indices.begin(); cur_it<inner_indices.end(); cur_it++)
    {

        //std::cout << "Cur inner index " << *cur_it << "column of lattice " << column(index(*t_begin)) << std::endl;
        *mat_columns=outer_val;
        //std::cout << "column set to " << *mat_columns << std::endl;
        mat_columns++;
        cur_row=row_map.begin();
        while (t_end-t_begin>0 && column(index_val(*t_begin))<*cur_it){

            t_begin++;
        }

        //less than greater than operators don't seem to work with tupleit (id3eally this function should just use index/data iterators instead.
        while (t_end-t_begin>0 && column(index_val(*t_begin))==*cur_it)
        {
            cur_row=std::lower_bound(cur_row,row_map.end(),row(index_val(*t_begin)));
            *mat_rows=cur_row-row_map.begin();
            //std::cout << "seting row to " << *mat_rows << std::endl;
            mat_rows++;
            outer_val++;
            t_begin++;


        }
        if(t_end-t_begin<1){

            break;
        }



    }

    while(mat_columns<=mat.outerIndexPtr()+mat.outerSize()){
        //std::cout << "final column " << mat_columns-mat.outerIndexPtr() << "set to " << mat.nonZeros() << std::endl;
        *mat_columns=mat.nonZeros();
        mat_columns++;
    }


}


template <class Derived>
inline void SparseLatticeBase<Derived>::to_matrix_rowmajor(const std::vector<index_type> &inner_indices,const std::vector<index_type> &column_map,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_rm& mat) const
{




    index_type outer_val=0;
    auto mat_columns=mat.innerIndexPtr();
    auto mat_rows=mat.outerIndexPtr();
    typename std::vector<index_type>::const_iterator cur_column;
    for(auto cur_it=inner_indices.begin(); cur_it<inner_indices.end(); cur_it++)
    {

        //std::cout << "Cur inner index " << *cur_it << "row of lattice " << row(index(*t_begin)) << std::endl;
        *mat_rows=outer_val;
        //std::cout << "row set to " << *mat_rows << std::endl;
        mat_rows++;
        cur_column=column_map.begin();


        while (t_end-t_begin>0 && row(index_val(*t_begin))<*cur_it)
            t_begin++;

        while (t_end-t_begin>0 && row(index_val(*t_begin))==*cur_it)
        {
            cur_column=std::lower_bound(cur_column,column_map.end(),column(index_val(*t_begin)));
            *mat_columns=cur_column-column_map.begin();
            //std::cout << "seting column to " << *mat_columns << std::endl;
            outer_val++;
            mat_columns++;
            t_begin++;
        }
        if(t_end-t_begin<1)
            break;


    }
    while(mat_rows<=mat.outerIndexPtr()+mat.outerSize()){
        //std::cout << "final row " << mat_rows-mat.outerIndexPtr() << "set to " << mat.nonZeros() << std::endl;
        *mat_rows=mat.nonZeros();
        mat_rows++;
    }


}

#ifdef LIBMIA_USE_SPARSE_SOLVE
template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(SparseLatticeBase<otherDerived> &b)
{

    this->check_solve_dims(b);



    typedef typename SparseSolveReturnType<Derived,otherDerived>::type c_type;
    typedef typename c_type::vector_type c_vector_type;

    typedef typename internal::index_type<otherDerived>::type b_index_type;
    typedef typename internal::storage_iterator<otherDerived>::type other_storage_iterator;
    typedef typename internal::index_iterator<otherDerived>::type other_index_iterator;


    //determine whether we want to binary search or just scan through elements
    bool a_search_flag=false, b_search_flag=false;
    if(this->depth()*std::log(this->size())<this->size())
        a_search_flag=true;
    if(this->depth()*std::log(b.size())<b.size())
        b_search_flag=true;

    otherDerived b_d=b.derived();
    Derived a_d=this->derived();
    sort(ColumnMajor); //tab/column major for A
    b.sort(ColumnMajor); //tab/column major for B

    //iterators for for indices and data
    index_iterator a_temp_begin=this->index_begin();
    auto a_index_end=this->index_end();
    index_iterator a_temp_end;
    auto b_temp_begin=b.index_begin();
    auto b_index_end=b.index_end();
    decltype(b_temp_begin) b_temp_end;


    //hold compressed indices for Sparse matrix of each tab
    std::vector<index_type> a_columns;
    //since all columns must have at least one non-zero entry to be fully-ranked, we can just do direct mapping to CSC
    a_columns.resize(this->width()+1);

    typedef typename Eigen::MappedSparseMatrix<data_type,Eigen::ColMajor,index_type> MappedSparseMatrix_cm; //Sparse Matrix type for A
    typedef Eigen::SuperLU<MappedSparseMatrix_cm> LU_decomp;    //LU decomposition type for A
    std::vector<typename MappedSparseMatrix_cm::Index> a_rows; //we make it size_t, incase index_type differs from size_t






    c_type c(this->width(),b.width(),this->depth());   //create dense lattice to return and allocate memory


    typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> mapped_vector; //set to *this's datatype - conversion happens if needed
    mapped_vector b_vector(b.height());




    for (int k=0; k<this->depth(); k++) //loop through every tab
    {


        //find the last occurence of current tab in both lattices
        a_temp_end=this->find_tab_end_idx(k,a_temp_begin,a_index_end,a_search_flag);
        //b could have an all-zero tab which is legal
        b_temp_begin=b.find_tab_start_idx(k,b_temp_begin,b_index_end,b_search_flag);
        //if all tabs in b are less than k, then we're done
        if(b_temp_begin==b_index_end)
            break;
        //if we can't find a tab equal to the current tab, we move on to the next candidate tab
        if(b.tab(*b_temp_begin)!=k)
            continue;

        b_temp_end=b.find_tab_end_idx(k,b_temp_begin,b_index_end,b_search_flag);

        //tab from a must have nonzeros
        if (a_temp_end!=a_temp_begin)
        {

            a_rows.resize(a_temp_end-a_temp_begin);
            //now we temporarily remap a_indices to rows
            a_rows[0]=this->row(*a_temp_begin); //set first to 0 row
            auto a_cur_it=a_temp_begin+1;
            auto a_rows_it=a_rows.begin()+1;
            size_t cur_column=0; //if size_t can't hold the number of columns, we're in big trouble anyway
            a_columns[0]=0;
            while(a_cur_it<a_temp_end){
                if(this->column(*a_cur_it)!=cur_column){
                    //check to make sure we didn't skip a row - if we did then we have a rank-deficient tab b/c of a zero-column
                    if(this->column(*a_cur_it)-cur_column>1){
                        std::stringstream t;
                        t << "Rank deficient tab. Column " << cur_column << " in tab " << k << " of LHS has zero entries.";
                        throw RankDeficientException(t.str());
                    }
                    cur_column++;
                    a_columns[cur_column]=a_cur_it-a_temp_begin;
                }
                *a_rows_it=this->row(*a_cur_it); //remap indices to rows
                a_cur_it++;
                a_rows_it++;
            }

            //end value of outer index array
            a_columns.back()=a_temp_end-a_temp_begin;

            //created a CCS matrix by mapping row and column vectors and also the pre-existing data of *this lattice
            MappedSparseMatrix_cm A=MappedSparseMatrix_cm(this->height(),this->width(),a_rows.size(),&a_columns[0],&a_rows[0],&(*(data_begin()+(a_temp_begin-this->index_begin())))); //map data to a compressed column matrix

            //compute LU decomposition
            LU_decomp lu_of_A(A);
            if(lu_of_A.info()!=Eigen::Success)
            {
                std::stringstream t;
                t << "Could not perform LU decomp on tab " << k << ".";
                throw RankDeficientException(t.str());
            }

            //loop through every column of b
            while(b_temp_begin<b_temp_end){
                b_vector.setZero(); //reset dense temp vector

                size_t b_cur_column=(size_t)b.column(*b_temp_begin);
                while(b_temp_begin<b_temp_end && b.column(*b_temp_begin)==b_cur_column){
                    b_vector(b.row(*b_temp_begin))=this->convert(b.derived().data_at(b_temp_begin-b.index_begin()));
                    b_temp_begin++;
                }


                //get wrapper for corresponding column of lattice c
                c_vector_type c_vector=c.column_vector(b_cur_column,k);
                //solve and store in lattice c
                lu_of_A._solve(b_vector,c_vector);
                if(lu_of_A.info()!=Eigen::Success)
                {
                    std::stringstream t;
                    t << "Solution process on tab " << k << " and column "<< b_cur_column << "of RHS failed.";
                    throw RankDeficientException(t.str());
                }


            }

        }
        else
        {
            std::stringstream t;
            t << "Rank deficient tab. Tab " << k << "has zero entries.";
            throw RankDeficientException(t.str());

        }
        //update current tab location
        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;

    }



    return c;
}
#else
template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(SparseLatticeBase<otherDerived> &b){
    //use delayed parsing, so the static assert will only trigger if the function is actually used within a compilation unit
    struct fake : std::false_type{};
    static_assert(fake::value,"You must define LIBMIA_USE_SPARSE_SOLVE (or turn it on if using CMake) and build and link to SuperLU to perform sparse solution of equations.");
}
#endif




/*! @} */




} //libMIA





#endif // SPARSELATTICEBASE_H
