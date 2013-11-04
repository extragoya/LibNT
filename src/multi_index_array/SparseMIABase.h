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


#include <boost/utility/enable_if.hpp>


#include "LibMIAException.h"
#include "LibMIAUtil.h"
#include "IndexUtil.h"
#include "MIA.h"
#include "MappedSparseLattice.h"
#include "LibMIAAlgorithm.h"
#include "FunctionUtil.h"

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
struct data_type_ref<SparseMIABase<Derived> >: public data_type_ref<Derived> {};

template<class Derived>
struct const_data_type_ref<SparseMIABase<Derived> >: public const_data_type_ref<Derived> {};

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

template<class Derived>
struct FinalDerived<SparseMIABase<Derived> >:public FinalDerived<Derived>{};


}



//!  Base class for sparse MIA classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  sparse multi-index array types. Provides operations and functions common to all sparse
  multi-index arrays.

  It assumes that non-zero data and corresponding indices are kept in two seperate data containers (as controlled by the Derived class).
  Indices are calculated using a linearized index. The lexographical precedence of indices used to caculate the linear index is controlled by mLinIdxSequence.
  By default the sequence is {0,1,2,...mOrder-1}, meaning calculation is idx0+dim0*idx1+dim0*dim1*idx2+.... However, mLinIdxSequence can be any
  permutation of the indices.

  \tparam Derived   should only be a sparse multi-index array class, eg, SparseMIA
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
    typedef typename internal::FinalDerived<SparseMIABase>::type FinalDerived;
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

    FinalDerived& final_derived() {

        return derived().final_derived();
    }
    /** \returns a const reference to the derived object */
    const FinalDerived& final_derived() const {

        return derived().final_derived();
    }

    //!Creates empty MIA of provided dimensionality. Number of variadic parameters must equal mOrder and be integer-like types
    template<typename... Dims>
    SparseMIABase(Dims... dims): MIA<SparseMIABase<Derived > >(dims...),mIsSorted(true) {
        init_linIdx_sequence();
    }


    //!Creates empty MIA of zero dimensionality
    SparseMIABase(): MIA<SparseMIABase<Derived > >(),mIsSorted(true) {
        init_linIdx_sequence();
    }




    //!Creates empty MIA of provided dimensionality
    SparseMIABase(const std::array<index_type,mOrder> &_dims,bool _is_sorted=true): MIA<SparseMIABase<Derived > >(_dims),mIsSorted(_is_sorted) {
        init_linIdx_sequence();
    }

    //! Returns iterator to data and indices at given linear index. Should check if it equals storage_end() before dereferencing
    /*!
        Will search the non-zero values for one whose linear index equals idx, meaning complexity is nnz*log(nnz). Note it assumes
        idx is calculated in the same way as mLinIdxSequence

    */
    const_storage_iterator atIdx(index_type idx) const{


        return find_idx(idx);

    }

    //!  Assignment based on given order. Will call Derived's operator
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up

        \param[in] otherMIA the other MIA
        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {3,1,2} this->at(1,2,3)==otherMIA.at(2,3,1).
                                Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,class index_param_type>
    void assign(const MIA<otherDerived>& otherMIA,const std::array<index_param_type,mOrder>& index_order){
        derived().assign(otherMIA.derived(),index_order);
    }


    //! Returns iterator to data and indices at given linear index. Should check if it equals storage_end() before dereferencing
    /*!
        Will search the non-zero values for one whose linear index equals idx, meaning complexity is nnz*log(nnz). Note it assumes
        idx is calculated in the same way as mLinIdxSequence

    */
    storage_iterator atIdx(index_type idx) {

        return find_idx(idx);
    }

    //! Iterator to the start of non-zeros entries. Deferencing iterator will provide a std::pair of references to data and index entries
    const_storage_iterator storage_begin() const
    {
        return iterators::makeTupleIterator(derived().data_begin(),derived().index_begin());

    }

    //! Iterator to the end of non-zeros entries.
    const_storage_iterator storage_end() const
    {
        return iterators::makeTupleIterator(derived().data_end(),derived().index_end());
    }
    //! Iterator to the start of non-zeros entries. Deferencing iterator will provide a std::pair of references to data and index entries
    storage_iterator storage_begin()
    {
        return iterators::makeTupleIterator(derived().data_begin(),derived().index_begin());

    }

    //! Iterator to the end of non-zeros entries.
    storage_iterator storage_end()
    {
        return iterators::makeTupleIterator(derived().data_end(),derived().index_end());
    }

    //! Iterator to the start of non-zero data container
    data_iterator data_begin()
    {


        return derived().data_begin();
    }

    //! Iterator to the end of non-zero data container
    data_iterator data_end()
    {


        return derived().data_end();
    }

    //! Iterator to the start of non-zero data container
    const_data_iterator data_begin() const
    {


        return derived().data_begin();
    }

    //! Iterator to the end of non-zero data container
    const_data_iterator data_end() const
    {


        return derived().data_end();
    }

    //! Iterator to the start of non-zero index container
    index_iterator index_begin()
    {


        return derived().index_begin();
    }

    //! Iterator to the end of non-zero index container
    index_iterator index_end()
    {


        return derived().index_end();
    }

    //! Iterator to the start of non-zero index container
    const_index_iterator index_begin() const
    {


        return derived().index_begin();
    }

    //! Iterator to the end of non-zero index container
    const_index_iterator index_end() const
    {


        return derived().index_end();
    }



    //! True if data and index containers are sorted based on the index container
    bool is_sorted() const{
        return mIsSorted;
    }



    //! Sort non-zero containers based on the given order
    /*!
        Will update the current sort order.

        \param[in]  _linIdxSequence the lexographical precedence to use in the sort - for instance {3,1,2} would sort based on the 3rd,1st, then 2nd index and update mLinIdxSequence accordingly
        \param[in] _stable whether to use stable sort or not. Unless there's good reason to, do not use stable sort, as its much slower (b/c is uses tuples of iterators)

    */
    void sort(const std::array<size_t,mOrder> & _linIdxSequence,bool _stable=false)
    {
        change_linIdx_sequence(_linIdxSequence);
        sort(_stable);
    }

    //! Sort non-zero containers based on the current sort order
    /*!

        \param[in] _stable whether to use stable sort or not. Unless there's good reason to, do not use stable sort, as its much slower (b/c is uses tuples of iterators)

    */
    void sort(bool _stable=false){
        sort(std::less<index_type>(),_stable);

    }

    //! Sort non-zero containers based on the current sort order that allows one to specify the comparison operation
    /*!

        \param[in] comp a binary predicate that returns true if its first argument is less than its second
        \param[in] _stable whether to use stable sort or not. Unless there's good reason to, do not use stable sort, as its much slower (b/c is uses tuples of iterators)

    */
    template<typename Comp>
    void sort(Comp comp,bool _stable=false)
    {


        if(!boost::is_same<Comp,std::less<index_type>>::value || !mIsSorted) //if comp is not std::less, then mIsSorted is meaningless
        {
            if(_stable){
                //stable_sort doesn't work with permute_iterator, b/c permute_iterator violates some of the iterator requirements
                //so we use tupleit, which is slower.
                std::stable_sort(storage_begin(),storage_end(),[] (const full_tuple& left,const full_tuple& right)
                {
                    return std::get<1>(left)<std::get<1>(right);
                } );
            }
            else{
                internal::Introsort(this->index_begin(),this->index_end(),comp,
                                    internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));
            }
        }
        if(boost::is_same<Comp,std::less<index_type>>::value) //if comp is not std::less, then mIsSorted should be false
            mIsSorted=true;


    }

    //!don't use - just here for benchmark purposes - will probably disappear in later versions
    void old_sort()
    {
        if(!mIsSorted)
        {
            std::sort(storage_begin(),storage_end(),[] (const full_tuple& left,const full_tuple& right)
                {
                    return std::get<1>(left)<std::get<1>(right);
                } );
        }
    }

    //!Prints non-zero values and indices
    void print()
    {

        this->sort(mDefaultLinIdxSequence);
        print_array(this->dims(),"Dimensions");
        std::cout << "Nonzeros " << this->size() << std::endl;
        std::cout << "Dimens" << std::endl;
        std::cout << "Index\t Indices\t Data" << std::endl;
        for(auto it=this->storage_begin();it<this->storage_end();++it){
            std::cout << index_val(*it) << "\t ";
            print_array_on_line(this->ind2sub(index_val(*it)));
            std::cout << "\t" << data_val(*it) << std::endl;
        }
        std::cout << std::endl;

    }

    //! Changes the lexicographical precedence used to calculate the linear indices.
    /*!
        If different than mLinIdxSequence, all index values are recalculated based upon _linIdx_sequence, and sorted flag will be set to false.
    */
    template<class index_param_type>
    void change_linIdx_sequence(const std::array<index_param_type,mOrder> & _linIdx_sequence)
    {
        static_assert(internal::check_index_compatibility<size_t,index_param_type>::type::value,"Must use an array convertable to index_type");
        if(_linIdx_sequence==mLinIdxSequence) //do nothing if we're already at the desired linIdxSequence
            return;

        for(auto& it: derived().m_indices){
            //std::cout << "idx " << it << std::endl;
            //print_array(this->ind2sub_reorder(it,mLinIdxSequence),"ind2sub");
            //std::cout << "sub2ind " << this->sub2ind_reorder(this->ind2sub_reorder(it,mLinIdxSequence),_linIdxSequence) << std::endl;

            //get the full set of indices from the current linear index, and then recalculate the linear index based on _linIdx_sequence
            it=this->sub2ind_reorder(this->ind2sub_reorder(it,mLinIdxSequence),_linIdx_sequence);
        }
        set_linIdx_sequence(_linIdx_sequence);
        mIsSorted=false;
    }

    //! Changes the lexicographical precedence stored in mLinIdxSequence, but does not actually update the actual index values. Only use this function if you are certain of the consequences.
    template<class index_param_type>
    void set_linIdx_sequence(const std::array<index_param_type,mOrder> & _linIdx_sequence)
    {
        static_assert(internal::check_index_compatibility<size_t,index_param_type>::type::value,"Must use an array convertable to index_type");
        mLinIdxSequence=_linIdx_sequence;

    }

    //! Resets the lexographical predence of the linear indices to the default precedence, i.e., {0,1...mOrder-1}. Index values are updated as well.
    void reset_linIdx_sequence(){
        change_linIdx_sequence(mDefaultLinIdxSequence);
    }

    //!converts a linear index calculated using mLinIdxSequence to a linear index calculated using mDefaultLinIdxSequence
    index_type convert_to_default_linIdxSequence(const index_type idx) const
    {


        return this->sub2ind_reorder(this->ind2sub_reorder(idx,mLinIdxSequence),mDefaultLinIdxSequence);
    }

    //!converts a linear index calculated using mDefaultLinIdxSequence to a linear index calculated using mLinIdxSequence
    index_type convert_from_default_sort(const index_type idx) const
    {

        //print_array(this->ind2sub(idx,mLinIdxSequence),"ind2sub");
        //std::cout << "sub2ind " << this->sub2ind(this->ind2sub(idx,mLinIdxSequence),mDefaultLinIdxSequence) << std::endl;
        return this->sub2ind_reorder(this->ind2sub_reorder(idx,mDefaultLinIdxSequence),mLinIdxSequence);
    }



    //!  Sets MIA index data to uniformly distributed random values within the valid index range.
    /*!
        Does not check that no duplicates are formed.

    */
    void rand_indices(){
        if(!this->dimensionality()) //do nothing is dimensionality is zero
            return;
        using namespace boost::numeric;

        boost::uniform_real<> uni_dist(0,this->dimensionality()-1);
        boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(LibMIA_gen(), uni_dist);
        typedef converter<index_type,boost::uniform_real<>::result_type,conversion_traits<index_type,boost::uniform_real<>::result_type>,def_overflow_handler,RoundEven<boost::uniform_real<>::result_type>> to_mdata_type;
        for (auto i=derived().index_begin();i<derived().index_end();++i){
            *i=to_mdata_type::convert(uni());
        }
        mIsSorted=false;

    }



    void setSorted(bool isSorted){
        mIsSorted=isSorted;
    }

    //!Returns the current lexicographical precedence used to compute the linearized indices
    const std::array<size_t,mOrder> & linIdxSequence() const
    {
        return mLinIdxSequence;
    }


    //!Comparison operator to another Sparse MIA. Will account for differing lexicographical precedences.
    template<class otherDerived>
    bool operator==(SparseMIABase<otherDerived>& otherMIA);

    template<class otherDerived>
    bool fuzzy_equals(SparseMIABase<otherDerived>& otherMIA, data_type precision);



    template<class otherDerived>
    bool operator==(const DenseMIABase<otherDerived>& otherMIA);

    template<class otherDerived>
    bool fuzzy_equals(const DenseMIABase<otherDerived>& otherMIA, data_type precision);

    template<class otherDerived>
    bool operator!=(const DenseMIABase<otherDerived>& otherMIA)
    {

        return !(*this==otherMIA.derived());

    }

    template<class otherDerived>
    bool operator!=(SparseMIABase<otherDerived>& otherMIA)
    {

        return !(*this==otherMIA.derived());

    }

    size_t size() const
    {
        return derived().size();
    }


    //! Performs destructive add (+=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is added like: this->at(x,y,z)+=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder. Should also assert a compile failure of the Derived sparse MIA does not
        allow resizing, eg MappedSparseMIA.
    */
    template<class otherDerived,typename index_param_type>
    SparseMIABase & plus_equal(MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){

        std::plus<data_type> op;
        derived().merge(b.derived(),op,index_order);
        return *this;
    }


    //! Performs destructive subtract (-=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is subtracted like: this->at(x,y,z)-=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder. Should also assert a compile failure of the Derived sparse MIA does not
        allow resizing, eg MappedSparseMIA.
    */
    template<class otherDerived,typename index_param_type>
    SparseMIABase & minus_equal(MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){
        std::minus<data_type> op;
        derived().merge(b.derived(),op,index_order);
        return *this;
    }

    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type plus_(SparseMIABase<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){

        std::plus<data_type> op;
        return outside_merge(b.derived(),op,index_order);

    }

    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type minus_(SparseMIABase<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){

        std::minus<data_type> op;
        return outside_merge(b.derived(),op,index_order);

    }

    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type plus_(const DenseMIABase<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order)const{


        typename MIAMergeReturnType<Derived,otherDerived>::type c;
        c.assign(b,index_order);
        return c.plus_equal(*this,internal::createAscendingIndex<mOrder>());


    }

    template<class otherDerived,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type minus_(const DenseMIABase<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order)const{

        typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
        typedef typename internal::data_type<CType>::type c_data_type;

        //do the inverse subtraction operation;

        CType c;
        c.assign(b,index_order);
        for(auto it=c.data_begin();it<c.data_end();++it)
            *it=-1*(*it);

        auto lambda=[](const c_data_type & a, const c_data_type & b){
            return a+b;
        };



        c.merge(*this,lambda,internal::createAscendingIndex<mOrder>());

        return c;

    }

    //statically delegates to different method is Op is std::minus
    template<typename otherDerived, typename Op,typename index_param_type,
        typename boost::disable_if<
            boost::is_same<
                Op,
                std::minus<
                    typename internal::data_type<
                        typename MIAMergeReturnType<
                            Derived,
                            otherDerived
                        >::type
                    >::type
                >
            >,
            int
        >::type=0
    >
    typename MIAMergeReturnType<Derived,otherDerived>::type outside_merge(const DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order) const{
        typename MIAMergeReturnType<Derived,otherDerived>::type c;

        c.assign(b,index_order);
        c.merge(*this,op,internal::createAscendingIndex<mOrder>());
        return c;

    }

    template<typename otherDerived, typename Op,typename index_param_type,
        typename boost::enable_if<
            boost::is_same<
                Op,
                std::minus<
                    typename internal::data_type<
                        typename MIAMergeReturnType<
                            Derived,
                            otherDerived
                        >::type
                    >::type
                >
            >,
            int
        >::type=0
    >
    typename MIAMergeReturnType<Derived,otherDerived>::type outside_merge(const DenseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order) const{

        return this->minus_(b,index_order);

    }

    template<typename otherDerived, typename Op>
    typename MIAMergeReturnType<Derived,otherDerived>::type outside_merge(const DenseMIABase<otherDerived> &b,const Op& op) const{
        return this->outside_merge(b,op,internal::createAscendingIndex<mOrder>());

    }

    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom.
    /*!
        For SparseMIAs, this function calls toLatticeSort, which will modify the ordering and indexing of data (but not it's actual data content)
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    MappedSparseLattice<data_type> toLatticeExpression(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices,bool columnSortOrder=true)
    {
        return toLatticeSort(row_indices, column_indices, tab_indices,columnSortOrder);

    }

    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom when the MIA is a temp object. Calls the same function as toLatticeExpression.
    /*!
        For SparseMIAs, this function calls toLatticeSort, which will modify the ordering and indexing of data (but not it's actual data content)
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    MappedSparseLattice<data_type> toLatticeDiscard(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices,bool columnSortOrder=true)
    {
        return toLatticeSort(row_indices, column_indices, tab_indices,columnSortOrder);

    }

    //! Flattens the MIA to a Lattice by creating a copy of the data.
    /*!
        \param[in] row_indices indices to map to the lattice rows - will perserve ordering
        \param[in] column_indices indices to map to the lattice columns - will perserve ordering
        \param[in] tab_indices indices to map to the lattice tabs - will perserve ordering
        \return MappedSparseLattice class that is a wrapper over the data and indices

        Will likely modify the data ordering of the SparseMIA
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
    MappedSparseLattice<data_type> toLatticeSort(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices,bool columnSortOrder);

    MappedSparseLattice<data_type> toStraightLattice(size_t number_of_row_indices, size_t number_of_column_indices,bool columnMajor=false);

    template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type noLatticeMult(SparseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

    template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type noLatticeMult(const DenseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

    index_type& index_val(const full_tuple & a)
    {
        return std::get<1>(a);

    }

    const index_type& index_val(const const_full_tuple & a) const
    {
        return std::get<1>(a);

    }


    data_type& data_val(const full_tuple & a)
    {
        return std::get<0>(a);

    }
    const data_type& data_val(const const_full_tuple & a) const
    {
        return std::get<0>(a);

    }

    //!returns the data value at the internal DOK data container index
    data_type& data_at(size_t dok_index)
    {
        return *(this->data_begin()+dok_index);

    }
    //!returns the data value at the internal DOK data container index
    const data_type& data_at(size_t dok_index) const
    {
        return *(this->data_begin()+dok_index);

    }

    template <class BinaryPredicate>
    index_iterator find_start_idx(index_iterator start_it, index_iterator end_it, const index_type & idx,BinaryPredicate predicate, bool search_flag=false);

    template <class BinaryPredicate>
    index_iterator find_end_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, BinaryPredicate predicate, bool search_flag=false);

    void set_linIdxSequence(const std::array<size_t,mOrder> & _linIdxSequence){
        mLinIdxSequence=_linIdxSequence;
    }

    template<typename otherDerived, typename Op,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type outside_merge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order);

    template<typename otherDerived, typename Op>
    typename MIAMergeReturnType<Derived,otherDerived>::type outside_merge(SparseMIABase<otherDerived> &b,const Op& op);

    //! Function for performing contraction and attraction. Best not to call directly, instead use the MIA algebra unary operation functionality, e.g., a(i,i)
    /*!
        \param[in] contract_indices list of indices undergoing a contraction
        \param[in] contract_partitions if more than one set of contractions is taking place, specifies how to partition contract_indices into corresponding sets of contractions
        \param[in] attract_indices list of indices undergoing an attraction
        \param[in] attract_partitions if more than one set of attractions is taking place, specifies how to partition attract_indices into corresponding sets of attractions
    */
    template<size_t no_con_indices,size_t no_con_partitions,size_t no_attract_indices,size_t no_attract_partitions>
    typename MIAUnaryType<Derived,internal::order<Derived>::value-no_con_indices-no_attract_indices+no_attract_partitions>::type contract_attract(const std::array<int,no_con_indices> & contract_indices,const std::array<int,no_con_partitions> & contract_partitions,const std::array<int,no_attract_indices> & attract_indices,const std::array<int,no_attract_partitions> & attract_partitions) ;

protected:

    //!keeps track of the order of dims used to calculate linear indices
    std::array<size_t,mOrder> mLinIdxSequence;
    //!keeps track of the order of dims used to calculate linear indices
    std::array<size_t,mOrder> mDefaultLinIdxSequence;
    //!keeps track of whether SparseMIA is sorted or not
    bool mIsSorted;

    storage_iterator find_idx(const index_type idx)
    {
        return std::lower_bound(storage_begin(),storage_end(),idx,[] (const full_tuple& _tuple,const index_type idx)
        {
            return std::get<1>(_tuple)<idx;
        } );

    }

    const_storage_iterator find_idx(const index_type idx)  const
    {
        return std::lower_bound(storage_begin(),storage_end(),idx,[] (const const_full_tuple& _tuple,const index_type idx)
        {
            return std::get<1>(_tuple)<idx;
        } );

    }

    template<class otherDerived, class Predicate>
    bool compare_with_dense(const DenseMIABase<otherDerived>& otherMIA,Predicate predicate);



    //!Comparison operator to another Sparse MIA. Will account for differing lexicographical precedences.
    template<class otherDerived, class Predicate>
    bool compare_with_sparse(SparseMIABase<otherDerived>& otherMIA,Predicate predicate)
    {
        if(this->dims()!=otherMIA.dims())
            return false;
        if(this->size()!=otherMIA.size())
            return false;
        this->sort();
        otherMIA.sort(mLinIdxSequence);
        auto it2=otherMIA.index_begin();
        for(auto it=this->index_begin();it<this->index_end();++it,++it2){
            if (*it!=*it2)
                return false;

        }
        auto data_it2=otherMIA.data_begin();
        for(auto data_it=this->data_begin();data_it<this->data_end();++data_it,++data_it2){
            if (!predicate(*data_it,*data_it2))
                return false;

        }
        return true;

    }

    void init_linIdx_sequence(){
        for(size_t i=0;i<mLinIdxSequence.size();++i)
            mLinIdxSequence[i]=i;
        mDefaultLinIdxSequence=mLinIdxSequence;
    }









    std::array<index_type,mOrder> ind2sub_reorder(index_type idx,const std::array<size_t,mOrder>& dim_order) const{
        //print_array(internal::ind2sub(idx, this->dims(),dim_order),"ind2sub");
        return internal::ind2sub_reorder(idx, this->dims(),dim_order);
    }

    index_type sub2ind_reorder(const std::array<index_type,mOrder> & indices,const std::array<size_t,mOrder>& dim_order) const{

        return internal::sub2ind_reorder(indices, dim_order,this->dims());
    }





    template<typename otherDerived, typename Op,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type outside_scanMerge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order);





private:

    friend class MIA<SparseMIABase>;

};



template<typename Derived>
template<typename otherDerived, typename Op>
typename MIAMergeReturnType<Derived,otherDerived>::type
SparseMIABase<Derived>::outside_merge(SparseMIABase<otherDerived> &b,const Op& op)
{

    return outside_merge(b,op,internal::createAscendingIndex<mOrder>());

}



template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type>
typename MIAMergeReturnType<Derived,otherDerived>::type
SparseMIABase<Derived>::outside_merge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order)
{
    this->check_merge_dims(b,index_order);
    std::array<size_t,mOrder> converted_index_order=array_converter<size_t>::convert(index_order);
    typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;


    if(b.is_sorted() || this->is_sorted()){
        //std::cout << "Performing Scan merge " <<std::endl;

        return outside_scanMerge(b,op,converted_index_order);
    }
    else{

        //std::cout << "Performing Sort merge " <<std::endl;
        CType C(*this);
        C.sortMerge(b,op,converted_index_order);
        return C;
    }

}

template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type>
typename MIAMergeReturnType<Derived,otherDerived>::type
SparseMIABase<Derived>::outside_scanMerge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order)
{
    typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
    if (b.is_sorted() && !this->is_sorted()){
        //get the order of lhs indices in terms of rhs
        auto lhsOrder=internal::reverseOrder(index_order);
        auto temp_linIdxSequence=this->mLinIdxSequence;
        internal::reorder_from(lhsOrder,b.linIdxSequence(),temp_linIdxSequence);
        this->sort(temp_linIdxSequence);
    }
    else if(this->is_sorted()&&!b.is_sorted()){
        auto b_linIdxSequence=b.linIdxSequence();
        internal::reorder_from(index_order,this->mLinIdxSequence,b_linIdxSequence);
        b.sort(b_linIdxSequence); //change b's sort order to it matches the index order, and also *this's current sort order

    }
    else if(this->is_sorted()&& b.is_sorted()){
        if (b.size()<this->size()){
            auto b_linIdxSequence=b.linIdxSequence();
            internal::reorder_from(index_order,this->mLinIdxSequence,b_linIdxSequence);
            b.sort(b_linIdxSequence); //change b's sort order to it matches the index order, and also *this's current sort order
        }
        else{
            //get the order of lhs indices in terms of rhs
            auto lhsOrder=internal::reverseOrder(index_order);
            auto temp_linIdxSequence=this->mLinIdxSequence;
            internal::reorder_from(lhsOrder,b.linIdxSequence(),temp_linIdxSequence);
            this->sort(temp_linIdxSequence);
        }
    }
    else
        throw MIAParameterException("Scan Merge should never have been called if both MIAs are unsorted");

    CType C(this->m_dims);
    C.change_linIdx_sequence(this->mLinIdxSequence);

    internal::outside_merge_sparse_storage_containers(C,*this,b,op);


    return C;
}

template<typename Derived>
auto SparseMIABase<Derived>::toStraightLattice(size_t number_of_row_indices, size_t number_of_column_indices,bool columnMajor) ->MappedSparseLattice<data_type>
{







    if(columnMajor){
        this->sort(mDefaultLinIdxSequence);
    }
    else{
        auto _linIdxSequence=mDefaultLinIdxSequence;
        for(size_t i=0;i<number_of_column_indices;++i){
            _linIdxSequence[i]=number_of_row_indices+i;
        }
        for(size_t i=number_of_column_indices;i<number_of_column_indices+number_of_row_indices;++i){
            _linIdxSequence[i]=i-number_of_column_indices;
        }
        sort(_linIdxSequence);
        //for now lattices don't actually change index values when doing column vs. row major - they just change the sort criteria, so we change it back to column major
        //but sorted row major
        change_linIdx_sequence(mDefaultLinIdxSequence);
    }




    index_type row_size=1, column_size=1, tab_size=1;


    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    for(size_t i=0;i<number_of_row_indices;++i)
    {

        row_size*=this->m_dims[i];
    }



    for(size_t i=number_of_row_indices;i<number_of_row_indices+number_of_column_indices;++i)
    {

        column_size*=this->m_dims[i];

    }


    for(size_t i=number_of_row_indices+number_of_column_indices;i<mOrder;++i)
    {

        tab_size*=this->m_dims[i];

    }

    //MappedSparseLattice<data_type> test=MappedSparseLattice<data_type>(derived().raw_data_ptr(),derived().raw_index_ptr(),this->size(),row_size,column_size,tab_size,columnMajor);
    return MappedSparseLattice<data_type>(derived().raw_data_ptr(),derived().raw_index_ptr(),this->size(),row_size,column_size,tab_size,columnMajor);






}


template<typename Derived>
template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
auto SparseMIABase<Derived>::toLatticeSort(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices,bool columnMajor) ->MappedSparseLattice<data_type>
{

    static_assert(internal::check_index_compatibility<index_type,idx_typeR>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeC>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeT>::type::value,"Must use an array convertable to index_type");
    std::array<size_t,mOrder> _linIdxSequence;

    if(columnMajor){
        _linIdxSequence=internal::concat_index_arrays(row_indices,column_indices,tab_indices);
        sort(_linIdxSequence);
    }
    else{
        _linIdxSequence=internal::concat_index_arrays(column_indices,row_indices,tab_indices);
        sort(_linIdxSequence);
        //for now lattices don't actually change index values when doing column vs. row major - they just change the sort criteria, so we change it back to column major
        //but sorted row major
        change_linIdx_sequence(internal::concat_index_arrays(row_indices,column_indices,tab_indices));
    }



    //statically check number of indices match up
    index_type row_size=1, column_size=1, tab_size=1;


    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    for(auto _row: row_indices)
    {

        row_size*=this->m_dims[_row];
    }



    for(auto _column: column_indices)
    {

        column_size*=this->m_dims[_column];

    }


    for(auto _tab: tab_indices)
    {

        tab_size*=this->m_dims[_tab];

    }

    //MappedSparseLattice<data_type> test=MappedSparseLattice<data_type>(derived().raw_data_ptr(),derived().raw_index_ptr(),this->size(),row_size,column_size,tab_size,columnMajor);
    return MappedSparseLattice<data_type>(derived().raw_data_ptr(),derived().raw_index_ptr(),this->size(),row_size,column_size,tab_size,columnMajor);


}

template<class Derived>
template<class otherDerived>
bool SparseMIABase<Derived>::operator==(const DenseMIABase<otherDerived>& otherMIA)
{
    typedef typename DenseMIABase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[](data_type a,other_data_type b){
        return a==b;
    };
    return compare_with_dense(otherMIA,pred);
}

template<class Derived>
template<class otherDerived>
bool SparseMIABase<Derived>::fuzzy_equals(const DenseMIABase<otherDerived>& otherMIA,data_type precision)
{

    typedef typename DenseMIABase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[precision](data_type a,other_data_type b){
        return isEqualFuzzy(a,b,precision);
    };
    return compare_with_dense(otherMIA,pred);
}


template<class Derived>
template<class otherDerived>
bool SparseMIABase<Derived>::operator==(SparseMIABase<otherDerived>& otherMIA)
{
    typedef typename SparseMIABase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[](data_type a,other_data_type b){
        return a==b;
    };
    return compare_with_sparse(otherMIA,pred);
}

template<class Derived>
template<class otherDerived>
bool SparseMIABase<Derived>::fuzzy_equals(SparseMIABase<otherDerived>& otherMIA,data_type precision)
{

    typedef typename SparseMIABase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[precision](data_type a,other_data_type b){
        return isEqualFuzzy(a,b,precision);
    };
    return compare_with_sparse(otherMIA,pred);
}

template<class Derived>
template<class otherDerived, class Predicate>
bool SparseMIABase<Derived>::compare_with_dense(const DenseMIABase<otherDerived>& otherMIA,Predicate predicate)
{

    //std::cout << "Entered == with dense" << std::endl;
    if(this->dims()!=otherMIA.dims()){
//        std::cout << "Dimension mismatch " << std::endl;
//        print_array(this->dims(),"this->dims()");
//        print_array(otherMIA.dims(),"otherMIA.dims()");
        return false;
    }
    this->sort();
    if (!this->size())
    {
        for(auto it=otherMIA.data_begin(); it<otherMIA.data_end(); ++it)
             if(!predicate(*it,0)){
                //std::cout << "Trigered not-zero no size " << it-otherMIA.data_begin() << " " << *it << std::endl;
                return false;
             }
        return true;
    }
    else
    {
        auto it=this->storage_begin();
        if (!predicate(otherMIA.atIdx(convert_to_default_linIdxSequence(index_val(*it))),data_val(*it))){
            //std::cout << "Trigered " << index_val(*it) << " " << convert_to_default_linIdxSequence(index_val(*it)) << " " << data_val(*it) << " " << otherMIA.atIdx(convert_to_default_linIdxSequence(index_val(*it))) << " " << data_val(*it)-otherMIA.atIdx(convert_to_default_linIdxSequence(index_val(*it))) << std::endl;
            return false;
        }

        for(index_type idx=0; idx<index_val(*(it)); idx++)
            if (!predicate(otherMIA.atIdx(convert_to_default_linIdxSequence(idx)),0)){
                //std::cout << "Trigered not-zero " << idx << " " << otherMIA.atIdx(convert_to_default_linIdxSequence(idx)) << std::endl;
                return false;
            }


        for(it=this->storage_begin()+1; it<this->storage_end(); ++it)
        {
            if (!predicate(otherMIA.atIdx(convert_to_default_linIdxSequence(index_val(*it))),data_val(*it)))
            {
                //std::cout << "Trigered " << index_val(*it) << " " << convert_to_default_linIdxSequence(index_val(*it)) << " " << data_val(*it) << " " << otherMIA.atIdx(convert_to_default_linIdxSequence(index_val(*it))) << " " << data_val(*it)-otherMIA.atIdx(convert_to_default_linIdxSequence(index_val(*it))) << std::endl;

                return false;
            }
            for(auto idx=index_val(*(it-1))+1; idx<index_val(*(it)); idx++)
                if (!predicate(otherMIA.atIdx(convert_to_default_linIdxSequence(idx)),0)){
                    //std::cout << "Trigered not-zero " << idx << " " << otherMIA.atIdx(convert_to_default_linIdxSequence(idx)) << std::endl;
                    return false;
                }

        }

        for(index_type idx=*(this->index_end()-1)+1; idx<this->m_dimensionality; idx++)
            if (!predicate(otherMIA.atIdx(convert_to_default_linIdxSequence(idx)),0)){
                //std::cout << "Trigered not-zero " << idx << " " << otherMIA.atIdx(convert_to_default_linIdxSequence(idx)) << std::endl;
                return false;
            }

        return true;
    }

}

template<typename Derived>
template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto  SparseMIABase<Derived>::noLatticeMult(SparseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type
{


    static_assert(Inter+L_outer==mOrder,"Both index arrays must index all indices in MIA");
    static_assert(Inter+R_outer==internal::order<otherDerived>::value,"Both index arrays must index all indices in MIA");
    typedef typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type RetType;
    typedef typename RetType::index_type c_index_type;
    typedef typename internal::index_type<otherDerived>::type b_index_type;
    std::array<index_type,Inter> l_inter_dims;
    std::array<index_type, L_outer> l_outer_dims;
    std::array<b_index_type,Inter> r_inter_dims;
    std::array<b_index_type,R_outer> r_outer_dims;
    //get inter and outer dimensionality and the individual dimensions that make up that number;
    index_type l_inter_size=internal::reorder_from(this->dims(), l_inter_idx,l_inter_dims);
    b_index_type r_inter_size=internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
    index_type l_outer_size=internal::reorder_from(this->dims(), l_outer_idx,l_outer_dims);
    b_index_type r_outer_size=internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);
    if(l_inter_dims.size()!=r_inter_dims.size() || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");

    std::array<c_index_type,L_outer+R_outer+Inter> c_dims;
    internal::concat_arrays(l_outer_dims,r_outer_dims,l_inter_dims,c_dims);
    RetType C (c_dims);
    C.reserve(this->size()*b.size()/std::pow(2,l_inter_dims.size())); //make an estimate of the number of elements in C based on how many indices are used for elemwise products
    if(Inter){
        change_linIdx_sequence(internal::concat_index_arrays(l_inter_idx,l_outer_idx)); //change order of linIdx calculation to inter then outer
        b.change_linIdx_sequence(internal::concat_index_arrays(r_inter_idx,r_outer_idx));
        auto lhs_compare=[&l_inter_size,this](const index_type &idx1, const index_type &idx2){
            return idx1%l_inter_size<idx2%l_inter_size;
        };
        this->sort(lhs_compare);
        auto rhs_compare=[&r_inter_size](const index_type &idx1, const index_type &idx2){
            return idx1%r_inter_size<idx2%r_inter_size;
        };
        b.sort(rhs_compare);
    }

    auto a_start_pred=[&](const index_type& idx, const index_type &val){
                return idx%l_inter_size<val;
    };
    auto a_end_pred=[&](const index_type& val, const index_type &idx){
                return val<idx%l_inter_size;
    };
    auto b_start_pred=[&](const b_index_type& idx, const b_index_type &val){
                return idx%r_inter_size<val;
    };
    auto b_end_pred=[&](const b_index_type& val, const b_index_type &idx){
                return val<idx%r_inter_size;
    };

    index_type cur_elem_wise_idx, cur_a_elem_wise_idx;
    b_index_type cur_b_elem_wise_idx;

    auto a_idx_begin=this->index_begin();
    auto b_idx_begin=b.index_begin();
    auto a_idx_end=this->index_end();
    auto b_idx_end=b.index_end();
    decltype(a_idx_end) cur_a_end;
    decltype(b_idx_end) cur_b_end;
    auto cur_a_idx=this->index_begin();
    auto cur_b_idx=b.index_begin();

    //determine whether we want to binary search or just scan through elements
    bool a_search_flag=false, b_search_flag=false;
    if(l_inter_size*log2(this->size())<this->size())
        a_search_flag=true;
    if(r_inter_size*log2(b.size())<b.size())
        b_search_flag=true;
    C.setSorted(true);
    while(cur_a_idx<a_idx_end && cur_b_idx<b_idx_end)
    {
        //if we're at the same tab, no work
        cur_a_elem_wise_idx=*cur_a_idx%l_inter_size;
        cur_b_elem_wise_idx=*cur_b_idx%r_inter_size;
        if(cur_a_elem_wise_idx==cur_b_elem_wise_idx)
        {
            cur_elem_wise_idx=cur_a_elem_wise_idx;
        }
        else if (cur_a_elem_wise_idx<cur_b_elem_wise_idx)  //if a's inter idx is less than b's
        {
            cur_elem_wise_idx=cur_b_elem_wise_idx;

            cur_a_idx=this->find_start_idx(cur_a_idx,a_idx_end,cur_elem_wise_idx,a_start_pred,a_search_flag);

            if(cur_a_idx==a_idx_end) //no inter idx in A is greater than or equal to B's current idx - so we're finished the entire multiplication routine
                break;
            //couldn't find an inter idx in A equal to B's current inter idx, but we found one greater than it - so now we need to try to find a matching inter idx in b
            if(*cur_a_idx%l_inter_size!=cur_elem_wise_idx)
                continue;
        }
        else
        {
            cur_elem_wise_idx=cur_a_elem_wise_idx;

            cur_b_idx=b.find_start_idx(cur_b_idx,b_idx_end,cur_elem_wise_idx,b_start_pred,b_search_flag);

            if(cur_b_idx==b_idx_end) //no inter idx in B is greater than or equal to A's current idx - so we're finished the entire multiplication routine
                break;
            //couldn't find an inter idx in B equal to A's current inter idx, but we found one greater than it - so now we need to try to find a matching inter idx in A
            if(*cur_b_idx%r_inter_size!=cur_elem_wise_idx)
                continue;
        }
        cur_a_end=this->find_end_idx(cur_a_idx,a_idx_end,cur_elem_wise_idx,a_end_pred,a_search_flag);
        cur_b_end=b.find_end_idx(cur_b_idx,b_idx_end,cur_elem_wise_idx,b_end_pred,b_search_flag);

        auto c_elemwise_idx=cur_elem_wise_idx*l_outer_size*r_outer_size; //the convention is to have indices in this order: [l_outer r_outer inter]
        for(auto a_it=cur_a_idx;a_it<cur_a_end;++a_it){
            auto cur_a_data_it=this->data_begin()+(a_it-a_idx_begin);
            auto l_outer_idx=(*a_it)/l_inter_size;
            for(auto b_it=cur_b_idx;b_it<cur_b_end;++b_it){
                auto cur_b_data_it=b.data_begin()+(b_it-b_idx_begin);
                auto r_outer_idx=(*b_it)/r_inter_size;
                C.push_back((*cur_a_data_it)*(*cur_b_data_it),l_outer_idx+r_outer_idx*l_outer_size+c_elemwise_idx);
            }
        }

        cur_a_idx=cur_a_end;
        cur_b_idx=cur_b_end;

    }

    return C;

}


template<typename Derived>
template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto  SparseMIABase<Derived>::noLatticeMult(const DenseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type
{


    static_assert(Inter+L_outer==mOrder,"Both index arrays must index all indices in MIA");
    static_assert(Inter+R_outer==internal::order<otherDerived>::value,"Both index arrays must index all indices in MIA");
    typedef typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type RetType;
    typedef typename RetType::index_type c_index_type;
    typedef typename internal::index_type<otherDerived>::type b_index_type;
    std::array<index_type,Inter> l_inter_dims;
    std::array<index_type, L_outer> l_outer_dims;
    std::array<b_index_type,Inter> r_inter_dims;
    std::array<b_index_type,R_outer> r_outer_dims;
    //get inter and outer dimensionality and the individual dimensions that make up that number;
    index_type l_inter_size=internal::reorder_from(this->dims(), l_inter_idx,l_inter_dims);
    internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
    index_type l_outer_size=internal::reorder_from(this->dims(), l_outer_idx,l_outer_dims);
    b_index_type r_outer_size=internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);
    if(l_inter_dims.size()!=r_inter_dims.size() || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");

    std::array<c_index_type,L_outer+R_outer+Inter> c_dims;
    internal::concat_arrays(l_outer_dims,r_outer_dims,l_inter_dims,c_dims);
    RetType C (c_dims);
    C.reserve(this->size()*b.dimensionality()*0.8); //make an estimate of the number of elements in C based on how many indices are used for elemwise products
    if(Inter){
        change_linIdx_sequence(internal::concat_index_arrays(l_inter_idx,l_outer_idx)); //change order of linIdx calculation to inter then outer
        std::function<bool(const index_type &idx_1, const index_type &idx_2)> lhs_compare=[&l_inter_size,this](const index_type &idx1, const index_type &idx2){
            return idx1%l_inter_size<idx2%l_inter_size;
        };
        this->sort(lhs_compare);
    }


    index_type cur_a_elem_wise_idx;
    b_index_type cur_b_elem_wise_idx;

    auto a_idx_begin=this->index_begin();
    auto a_idx_end=this->index_end();
    auto cur_a_idx=this->index_begin();

    C.setSorted(true);
    while(cur_a_idx<a_idx_end)
    {

        cur_a_elem_wise_idx=*cur_a_idx%l_inter_size;
        auto c_elemwise_idx=cur_a_elem_wise_idx*l_outer_size*r_outer_size; //the convention is to have indices in this order: [l_outer r_outer inter]
        cur_b_elem_wise_idx=internal::sub2ind(internal::ind2sub(cur_a_elem_wise_idx,r_inter_dims),r_inter_idx,b.dims());
        auto l_outer_idx=(*cur_a_idx)/l_inter_size;
        auto cur_a_data_it=this->data_begin()+(cur_a_idx-a_idx_begin);
        if(*cur_a_data_it){

            for(b_index_type b_cur_outer=0;b_cur_outer<r_outer_size;++b_cur_outer)
            {
                auto cur_b_data_it=b.data_begin()+cur_b_elem_wise_idx+internal::sub2ind(internal::ind2sub(b_cur_outer,r_outer_dims),r_outer_idx,b.dims());
                if(*cur_b_data_it){
                    C.push_back((*cur_a_data_it)*(*cur_b_data_it),l_outer_idx+b_cur_outer*l_outer_size+c_elemwise_idx);
                }

            }
        }


        cur_a_idx++;

    }

    return C;

}


template<typename Derived>
template<size_t no_con_indices,size_t no_con_partition,size_t no_attract_indices, size_t no_attract_partition>
typename MIAUnaryType<Derived,internal::order<Derived>::value-no_con_indices-no_attract_indices+no_attract_partition>::type
SparseMIABase<Derived>::contract_attract(const std::array<int,no_con_indices> & contract_indices,const std::array<int,no_con_partition> & contract_partition,const std::array<int,no_attract_indices> & attract_indices, const std::array<int,no_attract_partition> & attract_partition)
{

    typedef typename MIAUnaryType<Derived,internal::order<Derived>::value-no_con_indices-no_attract_indices+no_attract_partition>::type retType;

    //pull indices not undergoing an attraction or contraction
    auto copy_contract=internal::concat_index_arrays(contract_indices,attract_indices);
    std::sort(copy_contract.begin(),copy_contract.end());
    auto other_indices=internal::get_remaining_indices<size_t,no_con_indices+no_attract_indices,mOrder>(copy_contract);

    //get their dimensionality
    std::array<index_type,other_indices.size()> otherDims;
    auto other_dimensionality=internal::reorder_from(this->dims(), other_indices,otherDims);

    //sort a using other_indices as major, followed by attract indices and contract indices
    auto sort_order=internal::concat_index_arrays(contract_indices,attract_indices,other_indices);
    this->sort(sort_order);


    //get the dimensionality of indices involved in contraction
    index_type contract_dimensionality=index_type(1);
    size_t cur_idx=0;
    for(size_t i=0;i<no_con_partition;++i)
    {

       contract_dimensionality=internal::manual_int_power(this->dim(contract_indices[cur_idx]),contract_partition[i])*contract_dimensionality;
       cur_idx+=contract_partition[i];

    }

    //get the dimensionality of indices involved in attraction
    std::array<index_type,no_attract_partition> attract_index_ranges;
    cur_idx=0;
    index_type attract_dimensionality=1;
    for(size_t i=0;i<no_attract_partition;++i)
    {
       attract_index_ranges[i]= this->dim(attract_indices[cur_idx]);
       attract_dimensionality=internal::manual_int_power(attract_index_ranges[i],attract_partition[i])*attract_dimensionality;
       cur_idx+=attract_partition[i];

    }
    //add the attraction index ranges to the otherDims to get the returning dimensionality
    auto retDims=internal::concat_index_arrays(otherDims,attract_index_ranges);
    retType ret(retDims);


    index_type cur_non_zero_idx;
    index_type attract_idx=0;
    auto cur_it=this->index_begin();
    data_type sum;
    bool add_to_sum;
    //loop through all nonzeros of original MIA
    while(cur_it<this->index_end())
    {

        //get individual index locations in order of contract, attract, other_indices
        auto indices=internal::ind2sub(*cur_it, this->dims(),sort_order);

        add_to_sum=true;

        std::array<index_type,no_attract_partition> cur_attract_indices;
        size_t attract_start=no_con_indices;
        //if an attraction is taking place, examine the attraction indices of the current nonzero, and see if it falls on the appropriate diagonal
        for(size_t part_idx=0;part_idx<no_attract_partition;part_idx++)
        {
            auto check_idx=indices[attract_start]; //all indices in the current attraction set must be check_idx
            for(size_t _idx=attract_start+1;_idx<attract_start+attract_partition[part_idx];++_idx)
            {
                if(indices[_idx]!=check_idx){ //if it doesn't match throw a flag and break the loop
                    add_to_sum=false;
                    break;
                }

            }
            if(!add_to_sum)
                break;
            cur_attract_indices[part_idx]=check_idx; //update the array of attraction indices
            attract_start+=attract_partition[part_idx]; //caclulate the starting index of the next attraction set (if any)
        }
        if(add_to_sum){ //if all attraction indices match, calculate the linear index of them in the returning MIA
            attract_idx=internal::sub2ind(cur_attract_indices,attract_index_ranges);
        }
        else{ //if the attraction indices don't match, we are discarding this value, so just continue on to the next non-zero
            cur_it++;
            continue;
        }

        cur_non_zero_idx=*cur_it;
        sum=0;
        //to perform contraction examine all nonzeros that share the same non-contracted indices. Since we sorted the original MIA, all such nonzeros should be next to each other
        //if no contraction is being performed, this loop will just iterate once
        while(cur_it<this->index_end()&&cur_non_zero_idx/(contract_dimensionality)==(*cur_it)/(contract_dimensionality)){
            //std::cout << "Looking at index " << *cur_it << " default: " << convert_to_default_linIdxSequence(*cur_it) << std::endl;
            indices=internal::ind2sub(*cur_it, this->dims(),sort_order); //update our array of expanded indices
            //print_array(indices,"Current indices");
            add_to_sum=true;
            size_t contract_start=0;
            //examine each set of contractions, and if the indices for each set of contractions match, then include the data to the running sum
            for(size_t part_idx=0;part_idx<no_con_partition;part_idx++)
            {
                auto check_idx=indices[contract_start]; //all indices in the current contraction set must match check_idx
                for(size_t con_idx=contract_start+1;con_idx<contract_start+contract_partition[part_idx];++con_idx)
                {
                    if(indices[con_idx]!=check_idx){ //if the index doesn't match, this data element will be discarded, and we should continue on to the next one (if any)
                        add_to_sum=false;
                        break;
                    }

                }
                if(!add_to_sum)
                    break;
                contract_start+=contract_partition[part_idx]; //calculate the starting index of the next set of contractions
            }
            if(add_to_sum){
                //std::cout << "Adding to sum " <<std::endl;
                sum+=*(this->data_begin()+(cur_it-this->index_begin())); //if the contraction indices of the current data element all match, add it to the current sum

            }
            cur_it++; //examine the next data element
        }
        if(sum){
            //std::cout << "Pushing sum " <<sum << " to index " << cur_non_zero_idx/(contract_dimensionality*attract_dimensionality)+attract_idx*other_dimensionality <<std::endl;

            //push back the running sum, and also the index value
            ret.push_back(sum,cur_non_zero_idx/(contract_dimensionality*attract_dimensionality)+attract_idx*other_dimensionality);
        }

    }

    return ret;
    //print_array(other_indices,"other_indices");
    //print_array(retDims,"retDims");




}

template <class Derived>
template <class BinaryPredicate>
auto SparseMIABase<Derived>::find_start_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, BinaryPredicate pred,bool search_flag)->index_iterator
{


    if (search_flag){
        start_it= std::lower_bound(start_it,end_it,idx,pred);
    }
    else{

        while(start_it<end_it && pred(*start_it,idx))
                start_it++;

    }
    return start_it;

}

template <class Derived>
template <class BinaryPredicate>
auto SparseMIABase<Derived>::find_end_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, BinaryPredicate pred,bool search_flag)->index_iterator
{


    if (search_flag){
        start_it= std::upper_bound(start_it,end_it,idx,pred);
    }
    else{

        while(start_it<end_it && !pred(idx,*start_it))
            start_it++;

    }
    return start_it;

}



/*! @} */

}

#endif // SPARSEMIABASE_H_INCLUDED
