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

#include <type_traits>

#include <boost/utility/enable_if.hpp>

#include <boost/random/uniform_int_distribution.hpp>

#include "LibMIAException.h"
#include "LibMIAUtil.h"
#include "IndexUtil.h"
#include "MIA.h"
#include "MappedSparseLattice.h"
#include "LibMIAAlgorithm.h"
#include "FunctionUtil.h"
#include "LibMIATimSort.h"
#include "LibMIARadix.h"
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

/*!Helper structure to ensure that anytime the order used to calculated linear indices is changed, the dimension accumulator (used for breaking a linear index into its
     componenents) is automatically updated
*/
template<typename _index_type,size_t order>
struct linIdxContainer{
    std::array<size_t,order> m_linIdx;
    std::array<size_t,order> m_accumulator;
    const std::array<_index_type,order> & m_dims;

    linIdxContainer(const std::array<_index_type,order> & _linIdx, const std::array<_index_type,order> & _dims): m_linIdx(_linIdx),m_dims(_dims){

        calculate_accumulator();

    }

    linIdxContainer(const std::array<_index_type,order> & _dims): m_dims(_dims){

        m_linIdx.fill(0);
        m_accumulator.fill(0);

    }

    void calculate_accumulator(){
        m_accumulator[0]=1;
        for(size_t i=1;i<order;++i)
            m_accumulator[i]=m_accumulator[i-1]*m_dims[m_linIdx[i-1]];

    }
    linIdxContainer& operator=(const std::array<size_t,order> & _linIdx){
        m_linIdx=_linIdx;
        calculate_accumulator();
        return *this;
    }

    bool operator==(const std::array<size_t,order> & _linIdx)const{
        return m_linIdx==_linIdx;
    }
     bool operator!=(const std::array<size_t,order> & _linIdx)const{
        return !(m_linIdx==_linIdx);
    }
    _index_type& operator[](int idx) {
        return m_linIdx[idx];
    }
    const std::array<size_t,order> &  linIdxSequence() const{
        return m_linIdx;
    }
    std::array<size_t,order> &  linIdxSequence() {
        return m_linIdx;
    }
    const std::array<size_t,order> &  accumulator() const{
        return m_accumulator;
    }
    operator const std::array<size_t,order>()const{
        return m_linIdx;
    }
};

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
    typedef typename internal::data_type_ref<SparseMIABase>::type data_type_ref;
    typedef typename internal::const_data_type_ref<SparseMIABase>::type const_data_type_ref;
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
    typedef typename std::make_unsigned<index_type>::type cast_type; //type to use for casting to unsigned, for the purposes of performing division (faster if unsigned)
    constexpr static size_t mOrder=internal::order<SparseMIABase>::value;
    typedef typename std::array<size_t,mOrder> linIdxType;
    typedef typename MIA<SparseMIABase>::fast_divisor fast_divisor;
    typedef typename MIA<SparseMIABase>::unsigned_index_type unsigned_index_type;
    typedef typename MIA<SparseMIABase>::accumulator_type accumulator_type;
    typedef typename MIA<SparseMIABase>::fast_accumulator_type fast_accumulator_type;
    typedef typename MIA<SparseMIABase>::multiplier_type multiplier_type;

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
    SparseMIABase(Dims... dims): MIA<SparseMIABase<Derived > >(dims...),mIsSorted(true),mLinIdx(this->dims()) {
        init_linIdx_sequence();
    }


    //!Creates empty MIA of zero dimensionality
    SparseMIABase(): MIA<SparseMIABase<Derived > >(),mIsSorted(true),mLinIdx(this->dims()) {
        init_linIdx_sequence();
    }




    //!Creates empty MIA of provided dimensionality
    SparseMIABase(const std::array<index_type,mOrder> &_dims,bool _is_sorted=true): MIA<SparseMIABase<Derived > >(_dims),mIsSorted(_is_sorted),mLinIdx(this->dims()) {
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




    //! A sort test, to test sorting random indices using the different sorting options
    void sort_test(const int option=0){
        std::vector<index_type> scratch1;
        std::vector<data_type> scratch2;

        if (option==1){

            gfx::timsort(this->index_begin(),this->index_end(),this->data_begin(),this->data_end(),scratch1,scratch2);

            return;
        }
        else if(option==0){
            sort();
            return;
        }
        else if (option==7){


            scratch1.resize(this->size());
            scratch2.resize(this->size());
            internal::PerformPatienceRun(this->index_begin(),this->index_end(),this->data_begin(),scratch1,scratch2);
        }
        else if(option==8){
            scratch1.resize(this->size());
            scratch2.resize(this->size());


            internal::RadixSortStraight(this->index_begin(),this->data_begin(),scratch1.begin(),scratch2.begin(),this->size(),this->dimensionality());
        }
        else if(option==9){
            internal::RadixSortInPlace_PowerOf2Radix_Unsigned(this->index_begin(), this->data_begin(),this->size(),this->dimensionality() );

        }
        else if(option==10){
            scratch1.resize(this->size());
            scratch2.resize(this->size());
            internal::RadixSort(this->index_begin(),this->index_end(),this->data_begin(),scratch1.begin(), scratch2.begin(),this->dimensionality());
        }

    }

    //! used to test different permutation algorithms -
    void test_sort(const std::array<size_t,mOrder> & _linIdxSequence,const int option=0)
    {
        auto oldLinIdxSequence=this->linIdxSequence();

        change_linIdx_sequence(_linIdxSequence);

        std::vector<index_type> scratch1;
        std::vector<data_type> scratch2;

        if (option==1){

            gfx::timsort(this->index_begin(),this->index_end(),this->data_begin(),this->data_end(),scratch1,scratch2);

            return;
        }
        else if(option==2){
            sort();
            return;
        }
        else if(option==4){
            scratch1.reserve(this->size());
            scratch2.reserve(this->size());
            internal::NaturalMergesort(this->index_begin(),this->index_end(),this->data_begin(),scratch1,scratch2,std::less<index_type>());


            return;
        }

        else if(option==11){
            internal::RadixSortInPlace_PowerOf2Radix_Unsigned(this->index_begin(), this->data_begin(),this->size(),this->dimensionality() );
            return;
        }
        else if(option==12){

            scratch1.resize(this->size());
            scratch2.resize(this->size());


            internal::RadixSortStraight<2048,11,3000>(this->index_begin(),this->data_begin(),scratch1.begin(),scratch2.begin(),this->size(),this->dimensionality());


            return;
        }
        else if (option==0){
            scratch1.reserve(this->size());
            scratch2.reserve(this->size());
        }

        //if we performing a sort that tries to subdive the task into independent sublists
        auto shuffleSequence=internal::getShuffleSequence(oldLinIdxSequence,_linIdxSequence); //get the shuffle sequence from old to new

        auto reverseShuffleSequence=internal::getShuffleSequence(_linIdxSequence,oldLinIdxSequence); //get the shuffle sequence from new to old

        std::array<index_type,mOrder> divisor_list;
        divisor_list[0]=1;
        for(auto idx=1;idx<mOrder;++idx){
            divisor_list[idx]=divisor_list[idx-1]*this->dim(this->linIdxSequence()[idx-1]);
        }
        bool oldCanSkip=true;
        //find out which indices we can skip
        for(size_t curOldSequenceIdx=1;curOldSequenceIdx<mOrder;++curOldSequenceIdx){
            bool canSkip=true;

            for(int check_idx=mOrder-1;check_idx>0;--check_idx){
                if (shuffleSequence[check_idx]==curOldSequenceIdx)
                    break;
                else if(shuffleSequence[check_idx]<curOldSequenceIdx){
                    canSkip=false;
                    break;
                }
            }

            if(canSkip && !oldCanSkip){
                //std::cout << "Previous indices couldn't skip: Sorting indices below " << curOldSequenceIdx << std::endl;
                auto start_it=this->index_begin();


                fast_divisor _fast_divisor(divisor_list[curOldSequenceIdx]);
                auto cur_num_divisor=divisor_list[curOldSequenceIdx];
                //dedtermine region where we can sort
                while(start_it < this->index_end()){
                    auto end_it=start_it+1;
                    auto curValue=((unsigned_index_type)*start_it)/_fast_divisor;
                    auto nextValue=(curValue+1)*cur_num_divisor;
                    while(end_it<this->index_end() && *end_it<nextValue)
                        end_it++;


                        //std::cout << "start_it " << *start_it << " end_it " << *end_it << std::endl;

                        switch(option){
                            case 0:

                                internal::NaturalMergesort(start_it,end_it,this->data_begin()+(start_it-this->index_begin()),scratch1, scratch2,std::less<index_type>());
                                break;
                            case 3:
                                gfx::timsort(start_it,end_it,this->data_begin()+(start_it-this->index_begin()),this->data_begin()+(end_it-this->index_begin()),scratch1,scratch2);
                                break;
                            case 5:
                                internal::Introsort(start_it,end_it,std::less<index_type>(), internal::DualSwapper<index_iterator,data_iterator>(start_it,this->data_begin()+(start_it-this->index_begin())));
                                break;



                            case 7:
                                scratch1.resize(end_it-start_it);
                                scratch2.resize(end_it-start_it);
                                internal::PerformPatienceRun(start_it,end_it,this->data_begin()+(start_it-this->index_begin()),scratch1,scratch2);
                                break;
                        }

                    //



                    start_it=end_it;
                }
            }
            //if we've reached the end of indices to consider, then we need to just perform a straight sort
            else if(!canSkip && curOldSequenceIdx==mOrder-1){
                    //std::cout << "Reached the end: Sorting indices below " << curOldSequenceIdx << std::endl;
                    if(option==3)
                        gfx::timsort(this->index_begin(),this->index_end(),this->data_begin(),this->data_end(),scratch1,scratch2);
                    else if(option==5){
                        internal::Introsort(this->index_begin(),this->index_end(),std::less<index_type>(),
                                internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));
                    }
                    else if(option==0){

                        internal::NaturalMergesort(this->index_begin(),this->index_end(),this->data_begin(),scratch1,scratch2,std::less<index_type>());
                    }
                    else if (option==7){

                            scratch1.resize(this->size());
                            scratch2.resize(this->size());
                            internal::PerformPatienceRun(this->index_begin(),this->index_end(),this->data_begin(),scratch1,scratch2);
                    }

            }
            else{
               //std::cout << "Skipping index " << curOldSequenceIdx << std::endl;
            }
            oldCanSkip=canSkip;



        }





        this->setSorted(true);
    }

    //! Sort non-zero containers based on the given order
    /*!
        Will update the current sort order.

        \param[in]  _linIdxSequence the lexographical precedence to use in the sort - for instance {3,1,2} would sort based on the 3rd,1st, then 2nd index and update mLinIdxSequence accordingly


    */
    void sort(const std::array<size_t,mOrder> & _linIdxSequence)
    {


//        print_array(this->linIdxSequence(),"this->linIdxSequence()");
//        print_array(_linIdxSequence,"_linIdxSequence");
        typedef typename std::make_unsigned<index_type>::type unsigned_Type;
        if(!mIsSorted){

            change_linIdx_sequence(_linIdxSequence);
            this->sort();
        }
        else if(mOrder!=1 && _linIdxSequence!=this->linIdxSequence()){ //if *this is already sorted, and we're just changing the linIdxSequence, then we can do a more efficient merge-type sort


            auto oldLinIdxSequence=this->linIdxSequence();
            change_linIdx_sequence(_linIdxSequence); //change the linear indices to the new lexicographical precedence
            //auxiliary buffers
            std::vector<index_type> scratch1;
            std::vector<data_type> scratch2;
            //just perform introsort if the size is too small
            if(this->size()<3000){
               internal::Introsort(this->index_begin(),this->index_end(),std::less<index_type>(),internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));
                return;
            }



            auto reverseShuffleSequence=internal::getShuffleSequence(oldLinIdxSequence,_linIdxSequence); //get the shuffle sequence from new to old
			std::vector<unsigned_Type> divisors;
			std::vector<unsigned_Type> max_sizes;
			auto shuffled_dims = this->dims();
			internal::reorder_from(this->dims(), _linIdxSequence, shuffled_dims);
			bool first_stage = internal::setupPermute(reverseShuffleSequence, shuffled_dims, divisors, max_sizes);



            //create RadixShuffle object
			internal::RadixShuffle<index_type, data_type, 2048, 11, 3000> radixShuffle(max_sizes, divisors, this->dimensionality(), first_stage);
            //permute the sparse data based on the stage information provided
            radixShuffle.permute(this->index_begin(),this->data_begin(),this->size());


            return;

        }

        this->setSorted(true);
    }




     /*!locates the end of index sequence whose idx'th index (based on linIdxSequence) is equal to search_idx
        \param[in] do_search if do_search is true, will possibly an interpolation search. However, a search assumes ALL indices from start_it
                             to end_it are in ascending order based on the idx'th index, if this constraint is not true, the routine may not
                             return the correct result.
    */
    index_iterator find_end_specific_idx(index_type search_idx, size_t idx, index_iterator start_it, index_iterator end_it, bool do_search=false) const
    {

        if(start_it==end_it)
            return start_it;
        if(this->get_idx(*(end_it-1),idx)==search_idx)
            return end_it;

        if(do_search && (this->get_idx(*(end_it-1),idx)-this->get_idx(*start_it,idx))*log2((unsigned)(end_it-start_it))<end_it-start_it){
            start_it= internal::InterpolationSearchUpperBound(start_it, end_it,search_idx,
                [this,idx](const index_type & lhs){
                    return this->get_idx(lhs,idx);
                });
        }
        else{
            while(start_it<end_it && this->get_idx(*start_it,idx)<=search_idx)
                start_it++;
        }
        return start_it;


    }





    //returns the idx'th index of the linIdx, where linIdx is calculated based on current linIdxSequence
    index_type get_idx(index_type linIdx, size_t idx) const
    {

        //linear indices may be calculated in any shuffle order, e.g., {1,0,2}, so we must take that into account when pulling the i'th index
        typedef typename std::make_unsigned<index_type>::type cast_type;
        auto quotient=static_cast<cast_type>(linIdx)/static_cast<cast_type>(mLinIdx.accumulator()[idx]);
        return quotient%(static_cast<cast_type>(this->dim(this->linIdxSequence()[idx])));

    }



    //! Sort non-zero containers based on the current sort order
    /*!

        \param[in] _stable whether to use stable sort or not. Unless there's good reason to, do not use stable sort, as its much slower (b/c is uses tuples of iterators)

    */
    void sort(){
        sort(std::less<index_type>());

    }

    //! Sort non-zero containers based on the current sort order that allows one to specify the comparison operation
    /*!

        \param[in] comp a binary predicate that returns true if its first argument is less than its second
        \param[in] _stable whether to use stable sort or not. Unless there's good reason to, do not use stable sort, as its much slower (b/c is uses tuples of iterators)

    */
    template<typename Comp>
    void sort(Comp comp)
    {


        if(!boost::is_same<Comp,std::less<index_type>>::value || !mIsSorted) //if comp is not std::less, then mIsSorted is meaningless
        {



                internal::Introsort(this->index_begin(),this->index_end(),comp,
                                    internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));

        }
        if(boost::is_same<Comp,std::less<index_type>>::value) //if comp is not std::less, then mIsSorted should be false
            mIsSorted=true;


    }

    //!Prints non-zero values and indices
    void print(bool do_sort=true)
    {

//        if(!do_sort)
//            this->reset_linIdx_sequence();
//        else
//            this->sort(mDefaultLinIdxSequence);
        print_array(this->dims(),"Dimensions");
        std::cout << "Nonzeros " << this->size() << std::endl;
        std::cout << "Dimens" << std::endl;

        print_array(this->linIdxSequence(),"Linear Index Sequence");
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

        if(_linIdx_sequence==linIdxSequence()) //do nothing if we're already at the desired linIdxSequence
            return;
        if(mOrder==1)
            return;
        //setup fast reshuffle (based on static for)
        auto reorder_Dims=this->dims();
        auto new_reorder_Dims=this->dims();
        //get the dimensions ordered based on current linIdx
        internal::reorder_from(this->dims(),this->linIdxSequence(),reorder_Dims);
        internal::reorder_from(this->dims(),_linIdx_sequence,new_reorder_Dims);
        //get the shuffle sequence that suffles the new linIdx into the current one
        auto index_order=internal::getShuffleSequence(_linIdx_sequence,this->linIdxSequence());
        //create some helper values for linIdx shuffle


        accumulator_type dim_accumulator;
        fast_accumulator_type fast_dim_accumulator;
        multiplier_type multiplier;
        internal::create_shuffle_needs(reorder_Dims,new_reorder_Dims,index_order,dim_accumulator,fast_dim_accumulator,multiplier);



        for(auto& it: derived().m_indices){
            //std::cout << "idx " << it << std::endl;
            //print_array(this->ind2sub_reorder(it,mLinIdxSequence),"ind2sub");
            //std::cout << "sub2ind " << this->sub2ind_reorder(this->ind2sub_reorder(it,mLinIdxSequence),_linIdxSequence) << std::endl;

            //get the full set of indices from the current linear index, and then recalculate the linear index based on _linIdx_sequence
            //auto it2=this->sub2ind_reorder(this->ind2sub_reorder(it,linIdxSequence()),_linIdx_sequence);

            it=internal::reShuffleLinearIndex(it,multiplier,fast_dim_accumulator,dim_accumulator);





        }

        setLinIdxSequence(_linIdx_sequence);
        mIsSorted=false;
    }

    //! Changes the lexicographical precedence stored in mLinIdxSequence, but does not actually update the actual index values. Only use this function if you are certain of the consequences.
    template<class index_param_type>
    void setLinIdxSequence(const std::array<index_param_type,mOrder> & _linIdx_sequence)
    {
        static_assert(internal::check_index_compatibility<size_t,index_param_type>::type::value,"Must use an array convertable to index_type");
        mLinIdx=_linIdx_sequence;

    }

    //! Resets the lexographical predence of the linear indices to the default precedence, i.e., {0,1...mOrder-1}. Index values are updated as well.
    void reset_linIdx_sequence(){
        if (mOrder>1)
            change_linIdx_sequence(mDefaultLinIdxSequence);
    }



    //!converts a linear index calculated using mLinIdxSequence to a linear index calculated using mDefaultLinIdxSequence
    index_type convert_to_default_linIdxSequence(const index_type idx) const
    {


        return this->sub2ind_reorder(this->ind2sub_reorder(idx,this->linIdxSequence()),mDefaultLinIdxSequence);
    }

    //!converts a linear index calculated using mDefaultLinIdxSequence to a linear index calculated using mLinIdxSequence
    index_type convert_from_default_sort(const index_type idx) const
    {

        //print_array(this->ind2sub(idx,mLinIdxSequence),"ind2sub");
        //std::cout << "sub2ind " << this->sub2ind(this->ind2sub(idx,mLinIdxSequence),mDefaultLinIdxSequence) << std::endl;
        return this->sub2ind_reorder(this->ind2sub_reorder(idx,mDefaultLinIdxSequence),this->linIdxSequence());
    }



    //!  Sets MIA index data to uniformly distributed random values within the valid index range.
    /*!
        Does not check that no duplicates are formed.

    */
    void rand_indices(){
        if(!this->dimensionality()) //do nothing is dimensionality is zero
            return;
        using namespace boost::numeric;

        boost::random::uniform_int_distribution<index_type> uni_dist(0,this->dimensionality()-1);

        for (auto i=derived().index_begin();i<derived().index_end();++i){
            *i=uni_dist(LibMIA_gen());
        }
        mIsSorted=false;

    }



    void setSorted(bool isSorted){
        mIsSorted=isSorted;
    }

    //!Returns the current lexicographical precedence used to compute the linearized indices
    const std::array<size_t,mOrder> & linIdxSequence() const
    {
        return mLinIdx.linIdxSequence();
    }
    //!Returns the current lexicographical precedence used to compute the linearized indices
    std::array<size_t,mOrder> & linIdxSequence()
    {
        return mLinIdx.linIdxSequence();
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
    MappedSparseLattice<data_type> toLatticeSort(const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices,bool columnSortOrder=true);

    MappedSparseLattice<data_type> toStraightLattice(size_t number_of_row_indices, size_t number_of_column_indices,bool columnMajor=true);




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

    //!returns the data value corresponding to the given index iterator
    data_type_ref data_at(index_iterator index_it)
    {
        return *(this->data_begin()+(index_it-derived().index_begin()));

    }
    //!returns the data value corresponding to the given index iterator
    const_data_type_ref data_at(const_index_iterator index_it) const
    {
       return *(this->data_begin()+(index_it-derived().index_begin()));

    }

    //!returns the linear index value corresponding to the given index iterator.
    /*!This may be different than simply dereferencing the index iterator, as the values in the internal index array may be computed using a
       any arbitrary ordering scheme.
    */
    index_type index_at(const_index_iterator index_it) const
    {
        return this->convert_to_default_linIdxSequence(*index_it);

    }

    //!returns the index values corresponding to the given index iterator.
    /*!This may be different than simply dereferencing the index iterator and calling ind2sub, as the values in the internal index array may be computed using a
       any arbitrary ordering scheme.
    */
    std::array<index_type,mOrder> indices_at(const_index_iterator index_it) const
    {
        return this->ind2sub(this->convert_to_default_linIdxSequence(*index_it));

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


    //!checks whether the data is under the zero tolerance
    bool below_tolerance(const data_type & _data) const
    {
        if(std::abs(_data)<=SparseTolerance<data_type>::tolerance)
            return true;
        return false;
    }

protected:

    template<typename otherDerived, typename Op,typename index_param_type>
    typename MIAMergeReturnType<Derived,otherDerived>::type perform_outside_merge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order);


    template <class ElemGetter>
    index_iterator find_start_idx(index_iterator start_it, index_iterator end_it, const index_type & idx,const ElemGetter & getter, bool search_flag=false);

    template <class ElemGetter>
    index_iterator find_end_idx(index_iterator start_it, index_iterator end_it, const index_type & idx,const ElemGetter & getter, bool search_flag=false);


    index_iterator find_start_idx(index_iterator start_it, index_iterator end_it, const index_type & idx,  bool search_flag=false);


    index_iterator find_end_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, bool search_flag=false);



    //!Multiples *this with other MIA where the multiplication only consists of inter and outer products. NOTE, if &b==this, will perform a copy of b before performing multiplication
    template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type noLatticeMult(SparseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

    //!Multiples *this with other MIA where the multiplication only consists of inter and outer products. NOTE, if &b==this, will perform a copy of b before performing multiplication
    template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type perform_noLatticeMult(SparseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

    template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
    typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type noLatticeMult(const DenseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx);

    //!keeps track of whether SparseMIA is sorted or not
    bool mIsSorted;
    //!keeps track of the order of dims used to calculate linear indices
    linIdxContainer<index_type,mOrder> mLinIdx;
    //!keeps track of the order of dims used to calculate linear indices
    std::array<size_t,mOrder> mDefaultLinIdxSequence;


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

        this->sort();
        otherMIA.sort(this->linIdxSequence());
        auto it2=otherMIA.index_begin();
        auto it1=this->index_begin();
        while(it1<this->index_end() && it2<otherMIA.index_end()){
            auto & data1=this->data_at(it1);
            auto & data2=otherMIA.data_at(it2);
            if(this->below_tolerance(data1))
                it1++;
            else if(otherMIA.below_tolerance(data2))
                it2++;
            else{
                if (*it1!=*it2){
                    //std::cout << "Triggered Index " << *it1 << " " << *it2 << std::endl;
                    return false;
                }
                if (!predicate(data1,data2)){
                    //std::cout << "Triggered data " << data1 << " " << data2 << std::endl;
                    return false;
                }

                it1++;
                it2++;

            }

        }
        //in case nnz of each MIA is different, both index iterators will not have reached the end
        //so check remaining nnz are below the zero tolerance
        while(it1<this->index_end()){
            if(!this->below_tolerance(this->data_at(it1++))){
                //std::cout << "not below tolerance it1 " << this->data_at(it1) << " index " << *it1 << std::endl;
                return false;
            }
        }
        while(it2<otherMIA.index_end()){
            if(!otherMIA.below_tolerance(otherMIA.data_at(it2++))){
                //std::cout << "not below tolerance it2 " << otherMIA.data_at(it2) << " index " << *it2 << std::endl;
                return false;
            }
        }



        return true;

    }

    void init_linIdx_sequence(){
        for(size_t i=0;i<mOrder;++i)
            mDefaultLinIdxSequence[i]=i;
        this->setLinIdxSequence(mDefaultLinIdxSequence);
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


    template <class D1,class D2,bool D3,size_t D4> friend class MIA_Atom;
    template <class E1> friend class SparseMIABase;
    template <class E1> friend class DenseMIABase;
    template <class F1> friend class MIA;

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

    if(static_cast<void*>(&b)==static_cast<void*>(this)){ //since the mult uses sort, if b is the same MIA as *this, we must make a copy
        typedef typename internal::data_type<otherDerived>::type b_data_type;
        SparseMIA<b_data_type,internal::order<otherDerived>::value> b_copy(b);
        return this->perform_outside_merge(b_copy,op,index_order);
    }
    else
        return this->perform_outside_merge(b,op,index_order);
}

template<typename Derived>
template<typename otherDerived, typename Op,typename index_param_type>
typename MIAMergeReturnType<Derived,otherDerived>::type
SparseMIABase<Derived>::perform_outside_merge(SparseMIABase<otherDerived> &b,const Op& op,const std::array<index_param_type,internal::order<SparseMIABase>::value>& index_order){

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
    //print_array(index_order, "index_order in scanMerge");
    typedef typename MIAMergeReturnType<Derived,otherDerived>::type CType;
    if (b.is_sorted() && !this->is_sorted()){
        //get the order of lhs indices in terms of rhs
        auto lhsOrder=internal::reverseOrder(index_order);
        auto temp_linIdxSequence=this->linIdxSequence();
        internal::reorder_from(lhsOrder,b.linIdxSequence(),temp_linIdxSequence);
        this->sort(temp_linIdxSequence);
    }
    else if(this->is_sorted()&&!b.is_sorted()){

        auto b_linIdxSequence=b.linIdxSequence();


        internal::reorder_from(index_order,this->linIdxSequence(),b_linIdxSequence);

        b.sort(b_linIdxSequence); //change b's sort order to it matches the index order, and also *this's current sort order


    }
    else if(this->is_sorted()&& b.is_sorted()){
        if (b.size()<this->size()){
            auto b_linIdxSequence=b.linIdxSequence();
            internal::reorder_from(index_order,this->linIdxSequence(),b_linIdxSequence);

            b.sort(b_linIdxSequence); //change b's sort order to it matches the index order, and also *this's current sort order
        }
        else{
            //get the order of lhs indices in terms of rhs
            auto lhsOrder=internal::reverseOrder(index_order);
            auto temp_linIdxSequence=this->linIdxSequence();
            internal::reorder_from(lhsOrder,b.linIdxSequence(),temp_linIdxSequence);

            this->sort(temp_linIdxSequence);
        }
    }
    else
        throw MIAParameterException("Scan Merge should never have been called if both MIAs are unsorted");

    CType C(this->dims());
    C.change_linIdx_sequence(this->linIdxSequence());

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

    }




    index_type row_size=1, column_size=1, tab_size=1;


    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    for(size_t i=0;i<number_of_row_indices;++i)
    {

        row_size*=this->dim(i);
    }



    for(size_t i=number_of_row_indices;i<number_of_row_indices+number_of_column_indices;++i)
    {

        column_size*=this->dim(i);

    }


    for(size_t i=number_of_row_indices+number_of_column_indices;i<mOrder;++i)
    {

        tab_size*=this->dim(i);

    }
    this->setSorted(false);
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

    }



    //statically check number of indices match up
    index_type row_size=1, column_size=1, tab_size=1;


    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    for(auto _row: row_indices)
    {

        row_size*=this->dim(_row);
    }



    for(auto _column: column_indices)
    {

        column_size*=this->dim(_column);

    }


    for(auto _tab: tab_indices)
    {

        tab_size*=this->dim(_tab);

    }
    this->setSorted(false); //we don't know what will be done to the lattice data, so to be careful, we make sure we set sorted to false for the MIA.
                            //However there is no recourse if the user changes the lattice sort_order to something different.

    return MappedSparseLattice<data_type>(derived().raw_data_ptr(),derived().raw_index_ptr(),this->size(),row_size,column_size,tab_size,columnMajor,true);


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

        for(index_type idx=*(this->index_end()-1)+1; idx<this->dimensionality(); idx++)
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
    if(static_cast<void*>(&b)==static_cast<void*>(this)){ //since the mult uses sort, if b is the same MIA as *this, we must make a copy
        typedef typename internal::data_type<otherDerived>::type b_data_type;
        SparseMIA<b_data_type,internal::order<otherDerived>::value> b_copy(b);
        return this->perform_noLatticeMult(b_copy,l_inter_idx,l_outer_idx,r_inter_idx,r_outer_idx);
    }
    else
        return this->perform_noLatticeMult(b,l_inter_idx,l_outer_idx,r_inter_idx,r_outer_idx);
}



template<typename Derived>
template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto  SparseMIABase<Derived>::perform_noLatticeMult(SparseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
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
    this->sort(internal::concat_index_arrays(l_outer_idx,l_inter_idx)); //sort based on inter products, then outer products
    b.sort(internal::concat_index_arrays(r_outer_idx,r_inter_idx)); //sort based on inter products, then outer products

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

    auto a_elem_getter=[l_outer_size](index_type idx){
        return idx/l_outer_size;
    };
    auto b_elem_getter=[r_outer_size](b_index_type idx){
        return idx/r_outer_size;
    };

    //loop through each set of nonzeros simultaneously
    while(cur_a_idx<a_idx_end && cur_b_idx<b_idx_end)
    {
        //if we're at the same tab, no work
        cur_a_elem_wise_idx=a_elem_getter(*cur_a_idx);
        cur_b_elem_wise_idx=b_elem_getter(*cur_b_idx);
        if(cur_a_elem_wise_idx==cur_b_elem_wise_idx)
        {
            cur_elem_wise_idx=cur_a_elem_wise_idx;
        }
        else if (cur_a_elem_wise_idx<cur_b_elem_wise_idx)  //if a's inter idx is less than b's
        {
            cur_elem_wise_idx=cur_b_elem_wise_idx;

            cur_a_idx=this->find_start_idx(cur_a_idx,a_idx_end,cur_elem_wise_idx,a_elem_getter,true);

            if(cur_a_idx==a_idx_end) //no inter idx in A is greater than or equal to B's current idx - so we're finished the entire multiplication routine
                break;
            //couldn't find an inter idx in A equal to B's current inter idx, but we found one greater than it - so now we need to try to find a matching inter idx in b
            if(a_elem_getter(*cur_a_idx)!=cur_elem_wise_idx)
                continue;
        }
        else
        {
            cur_elem_wise_idx=cur_a_elem_wise_idx;

            cur_b_idx=b.find_start_idx(cur_b_idx,b_idx_end,cur_elem_wise_idx,b_elem_getter,true);

            if(cur_b_idx==b_idx_end) //no inter idx in B is greater than or equal to A's current idx - so we're finished the entire multiplication routine
                break;
            //couldn't find an inter idx in B equal to A's current inter idx, but we found one greater than it - so now we need to try to find a matching inter idx in A
            if(b_elem_getter(*cur_b_idx)!=cur_elem_wise_idx)
                continue;
        }
        cur_a_end=this->find_end_idx(cur_a_idx,a_idx_end,cur_elem_wise_idx,a_elem_getter,true);
        cur_b_end=b.find_end_idx(cur_b_idx,b_idx_end,cur_elem_wise_idx,b_elem_getter,true);

        auto c_elemwise_idx=cur_elem_wise_idx*l_outer_size*r_outer_size; //the convention is to have indices in this order: [l_outer r_outer inter]

        for(auto b_it=cur_b_idx;b_it<cur_b_end;++b_it){
            auto cur_b_data_it=b.data_begin()+(b_it-b_idx_begin);
            auto & cur_b_data=*cur_b_data_it;
            auto r_outer_idx=(*b_it)%r_outer_size;
            for(auto a_it=cur_a_idx;a_it<cur_a_end;++a_it){
                auto cur_a_data_it=this->data_begin()+(a_it-a_idx_begin);
                auto & cur_a_data=*cur_a_data_it;
                auto l_outer_idx=(*a_it)%l_outer_size;
                C.push_back(cur_a_data*cur_b_data,l_outer_idx+r_outer_idx*l_outer_size+c_elemwise_idx);
            }
        }

        cur_a_idx=cur_a_end;
        cur_b_idx=cur_b_end;

    }
    C.setSorted(true);
    return C;

}


template<typename Derived>
template<typename otherDerived, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto  SparseMIABase<Derived>::noLatticeMult(const DenseMIABase<otherDerived> &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<Derived,otherDerived,L_outer+R_outer+Inter>::type
{



    //std::cout << "Entered sparse * dense MIA " << std::endl;
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
    C.reserve(b.dimensionality()); //make an estimate of the number of elements in C based on how many indices are used for elemwise products


    this->sort(internal::concat_index_arrays(l_outer_idx,l_inter_idx)); //sort based on inter products, then outer products

//    if(Inter){ //if we inter products, sort *this based on inter product indices, with the outer product indices being in arbitrary order
//        change_linIdx_sequence(internal::concat_index_arrays(l_inter_idx,l_outer_idx)); //change order of linIdx calculation to inter then outer
//        std::function<bool(const index_type &idx_1, const index_type &idx_2)> lhs_compare=[&l_inter_size,this](const index_type &idx1, const index_type &idx2){
//            return idx1%l_inter_size<idx2%l_inter_size;
//        };
//        this->sort(lhs_compare);
//    }


    index_type cur_a_elem_wise_idx;
    b_index_type cur_b_elem_wise_idx;

    auto a_idx_begin=this->index_begin();
    auto a_idx_end=this->index_end();
    auto cur_a_idx=this->index_begin();


    while(cur_a_idx<a_idx_end)
    {

        cur_a_elem_wise_idx=*cur_a_idx/l_outer_size;
        auto c_elemwise_idx=cur_a_elem_wise_idx*l_outer_size*r_outer_size; //the convention is to have indices in this order: [l_outer r_outer inter]
        cur_b_elem_wise_idx=internal::sub2ind(internal::ind2sub(cur_a_elem_wise_idx,r_inter_dims),r_inter_idx,b.dims());
        auto l_outer_idx=(*cur_a_idx)%l_outer_size;
        auto cur_a_data_it=this->data_begin()+(cur_a_idx-a_idx_begin);
        if(*cur_a_data_it){

            for(b_index_type b_cur_outer=0;b_cur_outer<r_outer_size;++b_cur_outer)
            {
                auto cur_b_data_it=b.data_begin()+cur_b_elem_wise_idx+internal::sub2ind(internal::ind2sub(b_cur_outer,r_outer_dims),r_outer_idx,b.dims());
                if(*cur_b_data_it){
                    C.push_back((*cur_a_data_it)*(*cur_b_data_it),b_cur_outer +l_outer_idx*r_outer_size+c_elemwise_idx);
                }

            }
        }


        cur_a_idx++;

    }
    C.setLinIdxSequence(internal::concat_index_arrays(internal::createAscendingIndex<R_outer>(L_outer),internal::createAscendingIndex<L_outer>(0),internal::createAscendingIndex<Inter>(L_outer+R_outer)));
    C.setSorted(true);

    //std::cout << "Finished sparse * dense MIA " << std::endl;
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
    std::array<index_type,other_indices.size()+no_attract_partition> retDims;
    internal::concat_arrays(otherDims,attract_index_ranges,retDims);
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
template <class ElemGetter>
auto SparseMIABase<Derived>::find_start_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, const ElemGetter & getter,bool search_flag)->index_iterator
{


    if (search_flag && (*(end_it-1)-*start_it)*log2((unsigned)(end_it-start_it))<end_it-start_it){
        start_it= internal::InterpolationSearchLowerBound(start_it,end_it,idx,getter);
    }
    else{

        while(start_it<end_it && getter(*start_it)<idx)
            start_it++;

    }
    return start_it;


}

template <class Derived>
auto SparseMIABase<Derived>::find_start_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, bool search_flag)->index_iterator
{


    if (search_flag && (*(end_it-1)-*start_it)*log2((unsigned)(end_it-start_it))<end_it-start_it){
        start_it= internal::InterpolationSearchLowerBound(start_it,end_it,idx);
    }
    else{

        while(start_it<end_it && *start_it<idx)
            start_it++;

    }
    return start_it;

}



template <class Derived>
template <class ElemGetter>
auto SparseMIABase<Derived>::find_end_idx(index_iterator start_it, index_iterator end_it, const index_type & idx,const ElemGetter & getter, bool search_flag)->index_iterator
{


    if (search_flag && (*(end_it-1)-*start_it)*log2((unsigned)(end_it-start_it))<end_it-start_it){
        start_it= internal::InterpolationSearchUpperBound(start_it,end_it,idx,getter);
    }
    else{

        while(start_it<end_it && getter(*start_it) <=idx)
            start_it++;

    }
    return start_it;

}

template <class Derived>
auto SparseMIABase<Derived>::find_end_idx(index_iterator start_it, index_iterator end_it, const index_type & idx, bool search_flag)->index_iterator
{



    if (search_flag && (*(end_it-1)-*start_it)*log2((unsigned)(end_it-start_it))<end_it-start_it){
        start_it= internal::InterpolationSearchUpperBound(start_it,end_it,idx);
    }
    else{

        while(start_it<end_it && *start_it <=idx)
            start_it++;

    }
    return start_it;


}


////!Poly algorithm that chooses which type of sparse multiplication to perform based on whether operands are sparse, hyper-sparse, row-sparse, or column-sparse
//template<typename Derived1, typename Derived2, typename array_type,size_t Inner,size_t Inter,size_t L_outer,size_t R_outer>
//auto  sparseMIAMultPolyAlg(SparseMIABase<Derived1> &a,SparseMIABase<Derived2> &b,const std::array<array_type,L_outer>&l_row_idx,const std::array<array_type,Inner>&l_col_idx,const std::array<array_type,Inter>&l_tab_idx,const std::array<array_type,Inner>&r_row_idx,
//                           const std::array<array_type,R_outer>&r_col_idx,const std::array<array_type,Inter>&r_tab_idx)
//->typename MIAProductReturnType<Derived1,Derived2,L_outer+R_outer+Inter>::type *
//{
//
//    static_assert(L_outer+Inner+Inter==internal::order<Derived1>::value,"All indices must be accounted for when multiplying MIAs");
//    static_assert(R_outer+Inner+Inter==internal::order<Derived2>::value,"All indices must be accounted for when multiplying MIAs");
//
//    typedef typename MIAProductReturnType<Derived1,Derived2,L_outer+R_outer+Inter>::type retType;
//    typedef typename internal::index_type<retType>::type retIndexType;
//    retIndexType * cMIA;
//
//    auto l_row_size=dimensionality_from(a.dims(),l_row_idx);
//    auto inner_size=dimensionality_from(a.dims(),l_col_idx);
//    auto inter_size=dimensionality_from(a.dims(),l_tab_idx);
//    auto r_col_size=dimensionality_from(b.dims(),r_col_idx);
//
//    std::array<retIndexType,L_outer+R_outer+Inter> cMIA_dims;
//    size_t curIdx=0;
//    internal::reorder_from(a.dims(),l_row_idx,cMIA_dims,curIdx);
//    internal::reorder_from(b.dims(),r_col_idx,cMIA_dims,curIdx);
//    internal::reorder_from(a.dims(),l_tab_idx,cMIA_dims,curIdx);
//
//    auto a_per_tab_size=static_cast<double>(a.size())/inter_size; //valid if nonzeros distributed uniformly
//    auto b_per_tab_size=static_cast<double>(b.size())/inter_size;
//
//    bool l_row_sparse=(a_per_tab_size/l_row_size < 1);
//    bool l_col_sparse=(a_per_tab_size/inner_size < 1);
//    bool r_row_sparse=(b_per_tab_size/inner_size < 1);
//    bool r_col_sparse=(b_per_tab_size/r_col_size < 1);
//
//    //create the linIdxSequence if the resulting lattice ends up being sorted in RowMajor Order
//    auto _first =createAscendingIndex<R_outer>(L_outer);
//    auto _second=createAscendingIndex<L_outer>(0);
//    auto _third=createAscendingIndex<Inter>(L_outer+R_outer);
//    auto rowMajorLinIdx=std::array<size_t,L_outer+R_outer+Inter>;
//    internal::concat_arrays(_first,_second,_third,rowMajorLinIdx);
//
//    bool doOuter=false;
//    bool colMajor=true;
//    bool doAccum=true;
//
//
//    //both are hypersparse, perform outer-product algorithm
//    if((l_col_sparse&& r_row_sparse) ){
//        doOuter=true;
//    }
//    else if(!l_row_sparse && !l_col_sparse){ //if A is not hypersparse, but B is hypersparse in some way we just perform CSC mult
//        doAccum=true;
//        if(!l_row_sparse&&!l_col_sparse&&!r_row_sparse&&!r_col_sparse){ //exception is if both operands are NOT hypersparse, then choose one with smallest combo of nonzeros and dimenions
//            if (a.size()<b.size()){ //if |A|+number of rows is bigger than |B| and number of columns, it should be cheaper to perform CSR multiplication
//                colMajor=false;
//            }
//            else{
//                colMajor=true;
//            }
//        else{
//            colMajor=true;
//        }
//
//    }
//    else if(!r_row_sparse && !r_col_sparse){ //if B is not hypersparse, but A is hypersparse in some way we just perform CSR mult
//        colMajor=false;
//        doAccum=true;
//
//    }
//    else if(l_row_sparse && !l_col_sparse){ //if A is row sparse, and B is hypersparse in some way, we perform CSC mult but with no sparse accumulator
//        doAccum=false;
//        if(r_col_sparse && !r_row_sparse){ //exception is if B is row sparse
//            if (a.size()<b.size()){
//                colMajor=false;
//            }
//            else{
//                colMajor=true;
//            }
//        }
//        else{
//            colMajor=true;
//        }
//    }
//    else if(r_col_sparse && !r_row_sparse){ //if B is column sparse, and A is hypersparse in some way, we perform CSR mult but with no sparse accumulator
//        colMajor=false;
//        doAccum=false;
//    }
//
//    //Now we actually perform the multiplication
//
//    if(doOuter){ //doing an outer product mult
//        auto aLat=a.toLatticeExpression(l_row_idx,l_col_idx,l_tab_idx,true);
//        auto bLat=b.toLatticeExpression(r_row_idx,r_col_idx,r_tab_idx,false);
//        auto cLat=aLat.outer_times(bLat);
//        cMIA=new retType(cMIA_dims,std::move(cLat));
//    }
//    else if(colMajor){
//        auto aLat=a.toLatticeExpression(l_row_idx,l_col_idx,l_tab_idx,true); //column major
//        auto bLat=b.toLatticeExpression(r_row_idx,r_col_idx,r_tab_idx,true); //colum major
//        if(doAccum){
//            auto cLat=aLat.csc_times<false>(bLat);
//            cMIA=new retType(cMIA_dims,std::move(cLat));
//        }
//        else{
//            auto cLat=aLat.csc_no_accum<false>(bLat);
//            cMIA=new retType(cMIA_dims,std::move(cLat));
//        }
//    }
//    else{ //rowMajor
//        auto aLat=a.toLatticeExpression(l_row_idx,l_col_idx,l_tab_idx,false); //row major
//        auto bLat=b.toLatticeExpression(r_row_idx,r_col_idx,r_tab_idx,false); //row major
//        aLat.inPlaceTranspose(); //transpose them - they now both are column major
//        bLat.inPlaceTranspose();
//        if(doAccum){
//            auto cLat=bLat.csc_times<false>(aLat);
//            cLat.inPlaceTranspose(); //flip row and column dimensions around
//            cMIA =new retType(cMIA_dims,std::move(cLat),rowMajorLinIdx);
//        }
//        else{
//            auto cLat=bLat.csc_no_accum<false>(aLat);
//            cLat.inPlaceTranspose(); //flip row and column dimensions around
//            cMIA =new retType(cMIA_dims,std::move(cLat),rowMajorLinIdx);
//        }
//
//    }
//    return cMIA;
//
//
//}


/*! @} */

}

#endif // SPARSEMIABASE_H_INCLUDED
