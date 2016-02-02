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
#include <chrono>
#include <vector>
#include <type_traits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <array>
#include <unordered_map>
#include <type_traits>
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
//#include <boost/timer/timer.hpp>



#include <Eigen/Sparse>
#include <Eigen/Core>
//#include <sparsehash/dense_hash_map>

#include "Lattice.h"
#include "LibMIAAlgorithm.h"
#include "LibMIAUtil.h"
#include "LibMIARadix.h"
//#include "tupleit.hh" //tuple that supports simultaneous sorting of index and data arrays, based on the former. Quite slow, so ideally best to leave it out

//#define LM_CSC_TIMES 1 //define to perform old compressed column lattice multiplication - should be off
//#define LM_COLUMN_SEARCH 1 //define to perform column search during mult_scatter operation (not as efficient)
namespace LibMIA
{

namespace internal
{

template<class Derived>
struct Data<SparseLatticeBase<Derived> >: public Data<Derived> {};

template<class Derived>
struct Indices<SparseLatticeBase<Derived> >: public Indices<Derived> {};

template<class Derived>
struct data_type<SparseLatticeBase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct data_type_ref<SparseLatticeBase<Derived> >: public data_type_ref<Derived> {};

template<class Derived>
struct const_data_type_ref<SparseLatticeBase<Derived> >: public const_data_type_ref<Derived> {};

template<class Derived>
struct index_type<SparseLatticeBase<Derived> >: public index_type<Derived> {};

template<class Derived>
struct index_iterator<SparseLatticeBase<Derived> >: public index_iterator<Derived> {};

template<class Derived>
struct const_index_iterator<SparseLatticeBase<Derived> >: public const_index_iterator<Derived> {};

template<class Derived>
struct data_iterator<SparseLatticeBase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct const_data_iterator<SparseLatticeBase<Derived> >: public const_data_iterator<Derived> {};

template<class Derived>
struct full_iterator_tuple<SparseLatticeBase<Derived> >: public full_iterator_tuple<Derived> {};

template<class Derived>
struct const_full_iterator_tuple<SparseLatticeBase<Derived> >: public const_full_iterator_tuple<Derived> {};

template<class Derived>
struct const_full_tuple<SparseLatticeBase<Derived> >: public const_full_tuple<Derived> {};

template<class Derived>
struct full_tuple<SparseLatticeBase<Derived> >: public full_tuple<Derived> {};

template<class Derived>
struct storage_iterator<SparseLatticeBase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<SparseLatticeBase<Derived> >: public const_storage_iterator<Derived> {};

}


//!Data structure, whose ColumnSparse template parameter allows one to choose between different ways to access starting and end positions of a given column
/*!When we do column-by-column multiplication, typically it's assumed we have a CSC data structure that allows us to access each column in constant time
    e.g., see any publication on the CSC sparse matrix data structure. However, if you have many all-zero columns, the CSC data structure is wasteful,
    c.f., Buluc and Gilbert's publications on their doubly compressed sparse column data structure (DCSC). This helper structure will perform either CSC or
    DCSC column lookup depending on the ColumnSparse template parameter.
*/
template<class Derived,bool ColumnSparse>
struct MultHelper{



    MultHelper(SparseLatticeBase<Derived> & lat):mLat(lat){}

    ~MultHelper(){
        mLat.clearDCSCVectors();
    }
    void initMultHelpers(typename SparseLatticeBase<Derived>::index_iterator tab_begin,typename SparseLatticeBase<Derived>::index_iterator tab_end){
        mFastDivisor=mLat.createDCSCVectors(tab_begin,tab_end);
    }

    template<typename b_index_type>
    void getColumnOffset(b_index_type b_row2, typename SparseLatticeBase<Derived>::index_type & a_column_begin,typename SparseLatticeBase<Derived>::index_type & a_column_end) const{
        //search the column map for the start of the column equalling b's row

        decltype(*(mLat.JC.begin())) b_row=b_row2;
        auto _base=static_cast<unsigned_index_type>(b_row)/mFastDivisor;
        auto chunk_start = mLat.JC.begin()+mLat.AUX[static_cast<size_t>(_base)];
        auto chunk_end = mLat.JC.begin()+mLat.AUX[static_cast<size_t>(_base+1)];

        auto jc_it=chunk_start;
        //auto jc_it2=chunk_end;

        for(;jc_it<chunk_end;++jc_it){
            if(*jc_it>=b_row)
                break;
        }

        //if we didn't find the column, but found one greater than it or reached the end, update the current index_iterator of A
        if(jc_it ==chunk_end || *jc_it!=b_row){
            a_column_begin=0;
            a_column_end=0;
        }
        else{
            auto pos=jc_it-mLat.JC.begin();
            a_column_begin=mLat.CP[pos];
            a_column_end=mLat.CP[pos+1];
        }
    }


    SparseLatticeBase<Derived> & mLat;
    typename SparseLatticeBase<Derived>::fast_divisor mFastDivisor;
    typedef typename SparseLatticeBase<Derived>::unsigned_index_type unsigned_index_type;

};

template<class Derived>
struct MultHelper<Derived,false>{

    MultHelper(SparseLatticeBase<Derived> & lat):mLat(lat){}

    ~MultHelper(){
        mLat.clearCP();
    }

    void initMultHelpers(typename SparseLatticeBase<Derived>::index_iterator tab_begin,typename SparseLatticeBase<Derived>::index_iterator tab_end){
        mLat.createCP(tab_begin,tab_end);
    }

    template<typename b_index_type>
    void getColumnOffset(b_index_type b_row, typename SparseLatticeBase<Derived>::index_type & a_column_begin,typename SparseLatticeBase<Derived>::index_type & a_column_end) const{
        a_column_begin=mLat.CP[b_row];
        a_column_end=mLat.CP[b_row+1];
    }


    SparseLatticeBase<Derived> & mLat;

};


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


    typedef typename internal::data_type<SparseLatticeBase>::type data_type;
    typedef typename internal::data_type_ref<SparseLatticeBase>::type data_type_ref;
    typedef typename internal::const_data_type_ref<SparseLatticeBase>::type const_data_type_ref;
    typedef typename LibMIA::internal::index_type<SparseLatticeBase>::type index_type;
    typedef typename LibMIA::internal::Data<SparseLatticeBase>::type Data;
    typedef typename LibMIA::internal::Indices<SparseLatticeBase>::type Indices;
    typedef typename LibMIA::internal::full_iterator_tuple<SparseLatticeBase>::type full_iterator_tuple;
    typedef typename LibMIA::internal::const_full_iterator_tuple<SparseLatticeBase>::type const_full_iterator_tuple;
    typedef typename LibMIA::internal::full_tuple<SparseLatticeBase>::type full_tuple;
    typedef typename LibMIA::internal::const_full_tuple<SparseLatticeBase>::type const_full_tuple;
    typedef typename LibMIA::internal::storage_iterator<SparseLatticeBase>::type storage_iterator;
    typedef typename LibMIA::internal::const_storage_iterator<SparseLatticeBase>::type const_storage_iterator;
    typedef typename LibMIA::internal::index_iterator<SparseLatticeBase>::type index_iterator;
    typedef typename LibMIA::internal::const_index_iterator<SparseLatticeBase>::type const_index_iterator;
    typedef typename LibMIA::internal::data_iterator<SparseLatticeBase>::type data_iterator;
    typedef typename LibMIA::internal::const_data_iterator<SparseLatticeBase>::type const_data_iterator;
    typedef typename Eigen::SparseMatrix<data_type> SparseMatrix_cm;
    typedef typename Eigen::SparseMatrix<data_type,Eigen::RowMajor> SparseMatrix_rm;
    typedef typename Eigen::MappedSparseMatrix<data_type,Eigen::ColMajor,index_type> MappedSparseMatrix_cm;
    typedef typename Eigen::MappedSparseMatrix<data_type,Eigen::RowMajor,index_type> MappedSparseMatrix_rm;
    typedef typename std::make_unsigned<index_type>::type cast_type; //type to use for casting to unsigned, for the purposes of performing division (faster if unsigned)

	typedef typename Lattice<SparseLatticeBase>::fast_divisor fast_divisor;
	typedef typename Lattice<SparseLatticeBase>::unsigned_index_type unsigned_index_type;
	typedef typename Lattice<SparseLatticeBase>::accumulator_type accumulator_type;
	typedef typename Lattice<SparseLatticeBase>::fast_accumulator_type fast_accumulator_type;
	typedef typename Lattice<SparseLatticeBase>::multiplier_type multiplier_type;

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

	//old mult routine, that didn't use hash tables
	template <class otherDerived>
	typename SparseProductReturnType<Derived, otherDerived>::type old_times(SparseLatticeBase<otherDerived> &b);

    template <class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type operator*(const DenseLatticeBase<otherDerived> &b);


    //only enable for operands that have the same index_type
    template <class otherDerived>
    typename SparseSolveReturnType<Derived,otherDerived>::type solve(SparseLatticeBase<otherDerived> &b);

    template <class otherDerived>
    typename SparseSolveReturnType<Derived,otherDerived>::type solve(const DenseLatticeBase<otherDerived> &b);

    //only enable for operands that have the same index_type
    template <class otherDerived>
    typename SparseSolveReturnType<Derived,otherDerived>::type lsqr_solve(SparseLatticeBase<otherDerived> &b);

    template <class otherDerived>
    typename SparseSolveReturnType<Derived,otherDerived>::type lsqr_solve(const DenseLatticeBase<otherDerived> &b);

    data_type operator()(index_type _row, index_type _column, index_type _tab) const;





//
//

//




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


    void print(){
        print(this->size());
    }
    void print(index_type end_print)
    {

        assert(end_print<=this->size());
        assert(end_print>=0);
        auto _end=this->begin()+end_print;


        //sort(m_linIdxSequence);
        std::cout << "Val\t" ;
        std::cout << "Row\t"<< "Column\t" << "Tab\t" << "Index Val\n";
        for (storage_iterator i=this->begin(); i<_end; i++)
        {

            full_tuple temp=*i;
            //T temp2=std::get<0>(temp);

            std::cout << boost::tuples::get<0>(temp) <<"\t" ;
            std::cout << row(index_val(temp)) << "\t"<< column(index_val(temp)) << "\t" << tab(index_val(temp)) <<"\t" << index_val(temp) << std::endl;;
        }
    }

    index_type size() const
    {
        return derived().size();

    }
    bool idx_less(index_type a, index_type b) const
    {


        return a<b;

    }




    void sort(bool linIdxSequence=ColumnMajor)
    {





		
        if(!is_sorted()) //if *this is not sorted, we need to do a straight-up sort
        {

            //std::cout << "normal sort " << std::endl;
//             boost::timer::cpu_timer ltensor;
            if(m_linIdxSequence!=linIdxSequence)
                this->change_sort_order(linIdxSequence);

//            std::cout << "change sort order time " << boost::timer::format(ltensor.elapsed()) << std::endl;
//            ltensor=boost::timer::cpu_timer();
            internal::Introsort(this->index_begin(),this->index_end(),std::less<index_type>(),internal::DualSwapper<index_iterator,data_iterator>(this->index_begin(),this->data_begin()));
             //std::cout << "sort time " << boost::timer::format(ltensor.elapsed()) << std::endl;





        }
        else if(m_linIdxSequence!=linIdxSequence && is_sorted()){ //is *this is already sorted, but in a different sort order, we can take advantage of the partial ordering
            //std::cout << "special sort " << std::endl;

            this->change_sort_order(linIdxSequence);
           // boost::timer::cpu_timer ltensor;



			std::array<int, 3> reverseShuffleSequence = { 1, 0, 2 };//unlike MIAs, this will always be this value (get the shuffle sequence from new to old)
			std::vector<unsigned_index_type> divisors;
			std::vector<unsigned_index_type> max_sizes;
			std::array<index_type, 3> shuffle_dims = {this->height(),this->width(),this->depth()};
			if (linIdxSequence == RowMajor){
				shuffle_dims[0] = this->width();
				shuffle_dims[1] = this->height();
			}

			bool first_stage = internal::setupPermute(reverseShuffleSequence, shuffle_dims, divisors, max_sizes);



			//create RadixShuffle object
			internal::RadixShuffle<index_type, data_type, 2048, 11, 3000> radixShuffle(max_sizes, divisors, this->dimensionality(), first_stage);
			//permute the sparse data based on the stage information provided
			radixShuffle.permute(this->index_begin(), this->data_begin(), this->size());






        }
        this->set_sorted(true);


    }


    //!  Sets Lattice index data to uniformly distributed random values.
    /*!
    May cause duplicates

    */
    void rand_indices(){
        using namespace boost::numeric;

        boost::uniform_real<> uni_dist(0,this->dimensionality()-1);
        boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > uni(LibMIA_gen(), uni_dist);
        typedef converter<index_type,boost::uniform_real<>::result_type,conversion_traits<index_type,boost::uniform_real<>::result_type>,def_overflow_handler,RoundEven<boost::uniform_real<>::result_type>> to_mdata_type;
        for (auto i=derived().index_begin();i<derived().index_end();++i){
            *i=to_mdata_type::convert(uni());
        }
        m_is_sorted=false;

    }




    inline index_type row(index_type lin_index) const
    {

        if(this->linIdxSequence()==ColumnMajor){
			auto quotient = static_cast<unsigned_index_type>(lin_index) / this->height_divisor();
			return lin_index - quotient*this->height();

        }
		else{
			lin_index -= this->tab(lin_index)*this->height()*this->width();
			return static_cast<unsigned_index_type>(lin_index) / this->width_divisor();
		}

    }

    inline index_type column(index_type lin_index) const
    {
        if(this->linIdxSequence()==ColumnMajor){
			lin_index -= this->tab(lin_index)*this->height()*this->width();
			return static_cast<unsigned_index_type>(lin_index) / this->height_divisor();
        }
		else{
			auto quotient = static_cast<unsigned_index_type>(lin_index) / this->width_divisor();
			return lin_index - quotient*this->width();
		}

    }

    inline index_type tab(index_type lin_index) const
    {


		return static_cast<unsigned_index_type>(lin_index) / this->tab_divisor();


    }

    void set_sorted(bool _sorted)
    {

        m_is_sorted=_sorted;
    }


    bool is_sorted() const
    {

        return m_is_sorted;
    }
    bool linIdxSequence() const
    {

        return m_linIdxSequence;
    }

    //!Just changes the internal variables storing linIdxSequence and sorted flag - but doesn't alter any of the data
    void set_linIdxSequence(bool _linIdxSequence,bool _is_sorted=true){
        m_is_sorted=_is_sorted;
        m_linIdxSequence=_linIdxSequence;
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
		return boost::get<1>(a);

    }

    const index_type& index_val(const_full_tuple a) const
    {
        return boost::get<1>(a);

    }

    data_type& data_val(full_tuple a)
    {
		return boost::get<0>(a);

    }
    const data_type& data_val(const_full_tuple a) const
    {
		return boost::get<0>(a);

    }

	storage_iterator begin(){
		return derived().begin();
	}

	const_storage_iterator begin() const{
		return derived().begin();
	}

	storage_iterator end(){
		return derived().end();
	}

	const_storage_iterator end() const{
		return derived().end();
	}

    const data_type& data_at(size_t nnz_index) const
    {
        return *(data_begin()+nnz_index);

    }

    data_type& data_at(size_t nnz_index)
    {
        return *(data_begin()+nnz_index);

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







    void transpose(){
        inPlaceTranspose();
    }
    void inPlaceTranspose();


    template<class otherDerived>
    bool operator==(const DenseLatticeBase<otherDerived>& otherLat);

    //non constant b/c a sort may be involved
    template<class otherDerived>
    bool operator==(SparseLatticeBase<otherDerived>& otherLat);



    template<class otherDerived>
    bool fuzzy_equals(const DenseLatticeBase<otherDerived>& otherLat,data_type precision);

    //non constant b/c a sort may be involved
    template<class otherDerived>
    bool fuzzy_equals(SparseLatticeBase<otherDerived>& otherLat,data_type precision);

    bool below_tolerance(const data_type & _data) const;



    //non constant b/c a sort may be involved
    template<class otherDerived>
    bool operator!=(const DenseLatticeBase<otherDerived>& otherLat){
        return !(*this==otherLat);
    }

//    //!  Utility function for mapping tabs to column-major matrices to use in lattice mulitplication
//    /*!
//    \param[in] row_map A sorted vector mapping compressed indices to full row indices of the current tab
//    \param[in] inner_indices A sorted vector mapping shared column indices to full column indices of the current tab
//    \param[in] t_begin t_end Storage iterator to the nonzeros of the current tab
//    \param[out] mat Resulting CSC matrix. The CSC matrix should already have allocated vectors for data and indices.
//    */
//
//    //void to_matrix(const std::vector<index_type> &row_map,const std::vector<index_type> &inner_indices,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_cm& mat) const;
//
//
//    //!  Utility function for mapping tabs to row-major matrices to use in lattice mulitplication. See to_matrix for more information.
//
//    //void to_matrix_rowmajor(const std::vector<index_type> &inner_indices,const std::vector<index_type> &column_map,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_rm& mat) const;

    index_iterator find_tab_start_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag);
    index_iterator find_tab_end_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag);
    index_iterator find_tab_start_idx(index_type _tab,index_iterator start_it, index_iterator end_it);
    index_iterator find_tab_end_idx(index_type _tab,index_iterator start_it, index_iterator end_it);

    void change_sort_order(bool _linIdxSequence){


        if(this->linIdxSequence()==_linIdxSequence)
            return;
        if(_linIdxSequence==ColumnMajor)
            changeToColumnMajor();
        else
            changeToRowMajor();

        this->set_linIdxSequence(_linIdxSequence,false);

    }


    template <class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type exp_times(SparseLatticeBase<otherDerived> &b);


    template <class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type outer_times(SparseLatticeBase<otherDerived> &b);

    template <bool ColSparse, class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type csc_times(SparseLatticeBase<otherDerived> &b);

    template <bool ColSparse, class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type csc_no_accum(SparseLatticeBase<otherDerived> &b);

#ifndef LIBMIA_TEST_CPS //for testing purposes it's useful to be able to access these members
protected:
#endif
    std::vector<index_type> CP,JC,AUX,IR;
    void createCP(index_iterator tab_begin, index_iterator tab_end);

	void clearCP();

	void createJCCPIR(index_iterator tab_begin, index_iterator tab_end);

	void clearDCSCVectors();

protected:

    template<class b_index_type, class b_data_type,class ret_index_type,class super_data_type,bool RowSparse>
	bool mult_scatter(index_type a_tab_offset, b_index_type b_row, b_index_type b_column, b_data_type beta, std::vector<b_index_type>& restrict_libmia row_marker, std::vector<super_data_type> & restrict_libmia data_collector, std::vector<ret_index_type> & restrict_libmia c_indices, const MultHelper<Derived, RowSparse> & restrict_libmia multHelper);

	template<class b_index_type, class b_data_type, class ret_index_type,class Hasher>
	void mult_scatter(const index_iterator a_begin, const b_index_type b_row, const b_index_type b_column, const b_data_type beta, const fast_divisor& chunk_size,Hasher & row_hash,std::vector<ret_index_type> & c_indices) const;

    template<class b_index_type, class b_data_type, class ret_index_type,class ret_data_type,bool RowSparse>
	void mult_scatter(const index_iterator a_begin, const b_index_type b_row, const b_index_type b_column, const b_data_type beta, std::vector<ret_index_type> & restrict_libmia c_indices, std::vector<ret_data_type> & restrict_libmia c_data, const MultHelper<Derived, RowSparse> & restrict_libmia multHelper) const;

    void changeToColumnMajor(){

		std::array<unsigned_index_type, 3> reorder_Dims = { (unsigned_index_type)this->width(),(unsigned_index_type) this->height(), (unsigned_index_type)this->depth() };
		std::array<unsigned_index_type, 3> new_reorder_Dims = { (unsigned_index_type)this->height(),(unsigned_index_type) this->width(), (unsigned_index_type)this->depth() };

		//get the shuffle sequence that suffles the new linIdx into the current one
		std::array<int, 3> index_order = { 1, 0, 2 };
		//create some helper values for linIdx shuffle


		accumulator_type dim_accumulator;
		fast_accumulator_type fast_dim_accumulator;
		multiplier_type multiplier;
		internal::create_shuffle_needs(reorder_Dims, new_reorder_Dims, index_order, dim_accumulator, fast_dim_accumulator, multiplier);


		for (auto it = this->index_begin(); it < this->index_end(); ++it){
			*it = internal::reShuffleLinearIndex(*it, multiplier, fast_dim_accumulator, dim_accumulator);
		}
    }

    void changeToRowMajor(){

		std::array<unsigned_index_type, 3> reorder_Dims = { (unsigned_index_type)this->height(),(unsigned_index_type)this->width(), (unsigned_index_type)this->depth() };
		std::array<unsigned_index_type, 3> new_reorder_Dims = { (unsigned_index_type)this->width(), (unsigned_index_type)this->height(), (unsigned_index_type)this->depth() };

		//get the shuffle sequence that suffles the new linIdx into the current one
		std::array<int, 3> index_order = { 1, 0, 2 };
		//create some helper values for linIdx shuffle


		accumulator_type dim_accumulator;
		fast_accumulator_type fast_dim_accumulator;
		multiplier_type multiplier;
		internal::create_shuffle_needs(reorder_Dims, new_reorder_Dims, index_order, dim_accumulator, fast_dim_accumulator, multiplier);


		for (auto it = this->index_begin(); it < this->index_end(); ++it){
			*it = internal::reShuffleLinearIndex(*it, multiplier, fast_dim_accumulator, dim_accumulator);
		}
    }

    template<class otherDerived, class BinaryPredicate>
    bool compare_with_dense(const DenseLatticeBase<otherDerived>& otherLat,BinaryPredicate predicate);

    template<class otherDerived, class BinaryPredicate>
    bool compare_with_sparse(SparseLatticeBase<otherDerived>& otherLat,BinaryPredicate predicate);

    index_type full2lin_index(index_type _row, index_type _column, index_type _tab) const;

    void sparse_init(bool _is_sorted, bool _linIdxSequence)
    {

        m_is_sorted=_is_sorted;
        m_linIdxSequence=_linIdxSequence;


    }

    //! Returns the firstIdx - if lattice is sorted columnmajor, will return row, otherwise column
    inline index_type firstIdx(index_type _idx) const{
        if(this->linIdxSequence()==ColumnMajor){
            return (static_cast<cast_type>(_idx))%(static_cast<cast_type>(this->height()));
        }
        else
            return (static_cast<cast_type>(_idx))%(static_cast<cast_type>(this->width()));
    }
    //! Returns the second index - if lattice is sorted columnmajor, will return column, otherwise row
    inline index_type secondIdx(index_type _idx) const{
        if(this->linIdxSequence()==ColumnMajor){

            return ((static_cast<cast_type>(_idx)))/(static_cast<cast_type>(this->height()))%(static_cast<cast_type>(this->width()));

        }
        else{
            return ((static_cast<cast_type>(_idx)))/(static_cast<cast_type>(this->width()))%(static_cast<cast_type>(this->height()));
        }
    }





    fast_divisor createDCSCVectors(index_iterator tab_begin, index_iterator tab_end);

    

    fast_divisor createAUX(index_iterator tab_begin, index_iterator tab_end);



    




   

    template <class otherDerived, bool LSQR>
    typename SparseSolveReturnType<Derived,otherDerived>::type perform_solve(SparseLatticeBase<otherDerived> &b);

    template <class otherDerived,bool LSQR>
    typename SparseSolveReturnType<Derived,otherDerived>::type perform_solve(const DenseLatticeBase<otherDerived> &b);

    bool m_is_sorted;
    bool m_linIdxSequence;

private:
    template <class _Derived,bool RowSparse> friend class MultHelper;

};

template<typename Derived>
inline typename SparseLatticeBase<Derived>::index_type
SparseLatticeBase<Derived>::full2lin_index(index_type _row, index_type _column, index_type _tab) const
{

    if(this->linIdxSequence()==ColumnMajor)
        return _row+_column*this->height()+_tab*this->height()*this->width();
    else
        return _column+_row*this->width()+_tab*this->height()*this->width();

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




//!checks whether the data is under the zero tolerance
template<typename Derived>
bool SparseLatticeBase<Derived>::below_tolerance(const data_type & _data) const
{
    if(std::abs(_data)<=SparseTolerance<data_type>::tolerance)
        return true;
    return false;
}

//!This just swaps height and width values, and swaps current linIdxSequence. If you want a lattice with same linIdxSequence then you'll need to call change_sort_order.
template<typename Derived>
void SparseLatticeBase<Derived>::inPlaceTranspose()
{

    this->set_linIdxSequence(!this->linIdxSequence()); //update now switch back to old linIdxSequence, making a transpose
	auto temp= this->height();
	this->set_height(this->width());
	this->set_width(temp);



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
bool SparseLatticeBase<Derived>::operator==(SparseLatticeBase<otherDerived>& otherLat)
{
    typedef typename SparseLatticeBase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[](data_type a,other_data_type b){
        return a==b;
    };
    return compare_with_sparse(otherLat,pred);
}

template<class Derived>
template<class otherDerived>
bool SparseLatticeBase<Derived>::fuzzy_equals(SparseLatticeBase<otherDerived>& otherLat,data_type precision)
{
    typedef typename SparseLatticeBase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[precision](data_type a,other_data_type b){
        return isEqualFuzzy(a,b,precision);
    };
    return compare_with_sparse(otherLat,pred);
}

template<class Derived>
template<class otherDerived, class BinaryPredicate>
bool SparseLatticeBase<Derived>::compare_with_sparse(SparseLatticeBase<otherDerived>& otherLat,BinaryPredicate predicate)
{
    if(this->dims()!=otherLat.dims())
        return false;

    this->sort();
    otherLat.sort(this->linIdxSequence());
    auto it2=otherLat.index_begin();
    auto it1=this->index_begin();
    while(it1<this->index_end() && it2<otherLat.index_end()){
        auto & data1=this->data_at(it1);
        auto & data2=otherLat.data_at(it2);
        if(this->below_tolerance(data1))
            it1++;
        else if(otherLat.below_tolerance(data2))
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
    while(it1++<this->index_end()){
        if(!this->below_tolerance(this->data_at(it1))){
            //std::cout << "not below tolerance it1 " << this->data_at(it1) << " index " << *it1 << std::endl;
            return false;
        }
    }
    while(it2++<otherLat.index_end()){
        if(!otherLat.below_tolerance(otherLat.data_at(it2))){
            //std::cout << "not below tolerance it2 " << otherMIA.data_at(it2) << " index " << *it2 << std::endl;
            return false;
        }
    }



    return true;


}

template<class Derived>
template<class otherDerived>
bool SparseLatticeBase<Derived>::fuzzy_equals(const DenseLatticeBase<otherDerived>& otherLat,data_type precision)
{
    typedef typename DenseLatticeBase<otherDerived>::data_type other_data_type;
    std::function<bool(data_type,other_data_type)> pred=[precision](data_type a,other_data_type b){
        return isEqualFuzzy(a,b,precision);
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
            if(!predicate(*it,0))
                return false;

        return true;
    }
    else{

        this->sort(ColumnMajor);
        auto it=this->begin();
        if (!predicate(otherLat.atIdx(index_val(*it)),data_val(*it))){
            //std::cout << "Trigered " << index_val(*it) << " " << data_val(*it) << " " << otherLat.atIdx(index_val(*it)) << " " << data_val(*it)-otherLat.atIdx(index_val(*it)) << std::endl;
            return false;
        }

        for(index_type idx=0; idx<index_val(*(it)); idx++)
             if(!predicate(otherLat.atIdx(idx),0)){
                //std::cout << "Trigered not-zero " << idx << " " << otherLat.atIdx(idx) << std::endl;
                return false;
             }


        for(it=this->begin()+1; it<this->end(); ++it)
        {
            if (!predicate(otherLat.atIdx(index_val(*it)),data_val(*it)))
            {
                //std::cout << "Trigered " << index_val(*it) << " " << data_val(*it) << " " << otherLat.atIdx(index_val(*it)) << " " << data_val(*it)-otherLat.atIdx(index_val(*it)) << std::endl;

                return false;
            }
            for(auto idx=index_val(*(it-1))+1; idx<index_val(*(it)); idx++)
                if(!predicate(otherLat.atIdx(idx),0)){
                    //std::cout << "Trigered not-zero " << idx << " " << otherLat.atIdx(idx) << std::endl;
                    return false;
                }

        }

        for(index_type idx=*(this->index_end()-1)+1; idx<this->height()*this->width()*this->depth(); idx++)
             if(!predicate(otherLat.atIdx(idx),0)){
                //std::cout << "Trigered not-zero " << idx << " " << otherLat.atIdx(idx) << std::endl;
                return false;
            }

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
auto SparseLatticeBase<Derived>::find_tab_start_idx(index_type _tab,index_iterator start_it, index_iterator end_it)->index_iterator
{
    if(start_it==end_it)
        return start_it;
    if(this->tab(*start_it)==_tab)
        return start_it;
    bool search_flag=false;
    if(log2((unsigned)(end_it-start_it))< (end_it-start_it)/(this->tab(*(end_it-1))-this->tab(*start_it)))
        search_flag=true;
    return find_tab_start_idx(_tab,start_it,end_it,search_flag);
//if(this->depth()*log2((unsigned)(end_it-start_it))<this->size())
}

template <class Derived>
auto SparseLatticeBase<Derived>::find_tab_end_idx(index_type _tab,index_iterator start_it, index_iterator end_it, bool search_flag)->index_iterator
{

    if(search_flag)
        start_it= internal::InterpolationSearchUpperBound(start_it, end_it,_tab,
            [this](const index_type & lhs){
                return this->tab(lhs);
            }
        );
//        std::upper_bound(start_it,end_it,_tab,[this](index_type lhs, index_type rhs){
//                return lhs<this->tab(rhs);
//            });
    else{
        while(start_it<end_it && this->tab(*start_it)<=_tab)
            start_it++;
    }
    return start_it;

}

template <class Derived>
auto SparseLatticeBase<Derived>::find_tab_end_idx(index_type _tab,index_iterator start_it, index_iterator end_it)->index_iterator
{
    if(start_it==end_it)
        return start_it;
    if(this->tab(*(end_it-1))==_tab)
        return end_it;
    bool search_flag=false;
    if(log2((unsigned)(end_it-start_it))< (end_it-start_it)/(this->tab(*(end_it-1))-this->tab(*start_it)))
        search_flag=true;
    return find_tab_end_idx(_tab,start_it,end_it,search_flag);


}



template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::operator*(const DenseLatticeBase<otherDerived> &b){

	//TODO, right now this function doesn't check whether *this is row-sparse. If so, then the code execute as below. However, if not, then a sparse accumulator can be used. Should check for this
    //std::cout << "Entered sparse*dense " << std::endl;
#ifdef LIBMIA_CHECK_DIMS
    this->check_mult_dims(b);
#endif

	typedef typename SparseProductReturnType<Derived, otherDerived>::type CType;
	typedef typename internal::data_type<CType>::type c_data_type;
	typedef typename internal::index_type<CType>::type c_index_type;
	CType ret(this->height(), b.width(), this->depth());
	ret.reserve(std::min(b.dimensionality(), this->size()));
    //std::cout << "Finished mult " << std::endl;

	if (!this->size()){
		return ret;
	}
	this->sort(ColumnMajor);
	
	c_data_type cur_c_data;
	c_index_type c_column_begin = 0;
	c_index_type c_column_size = 0;
	index_iterator a_cur_it;
	index_iterator a_tab_begin=this->index_begin();
	index_iterator a_index_end = this->index_end();
	index_iterator a_temp_end;
	index_type a_cur_column;
	size_t old_c_size;
	bool a_search_flag = false;
	if (this->depth()*log2(this->size())<this->size())
		a_search_flag = true;
	
	while (a_tab_begin < a_index_end){ //we iterate through each tab of A
		index_type cur_tab = this->tab(*a_tab_begin); //get the current tab
		a_temp_end = this->find_tab_end_idx(cur_tab, a_tab_begin, a_index_end, a_search_flag); //find the end of current tab
		for (index_type b_cur_column = 0; b_cur_column < b.width(); ++b_cur_column){ //now iterate through each column of B			
			a_cur_it = a_tab_begin;
			old_c_size = ret.size();
			while (a_cur_it < a_temp_end){ //for each column of B, we iterate completely through A, and accumulate values for corresponding column of C
				a_cur_column = this->column(*a_cur_it);
				
				while (a_cur_it < a_temp_end && this->column(*a_cur_it) == a_cur_column){ //collect non-zeros of C arising form cur column of A
					cur_c_data = this->data_at(a_cur_it)*b(a_cur_column, b_cur_column, cur_tab);
					if (std::abs(cur_c_data)>SparseTolerance<double>::tolerance){
						ret.push_back(cur_c_data, this->row(*a_cur_it));						
					}
					a_cur_it++;
				}
				
				
				
			}
			//now sort and sum duplicate non-zero values, producing the resulting column of C
			internal::RadixSortInPlace_PowerOf2Radix_Unsigned(ret.index_begin() + old_c_size, ret.data_begin() + old_c_size, ret.size() - old_c_size, this->height());
			auto diff = collect_duplicates_function(ret.index_begin() + old_c_size, ret.index_end(), ret.data_begin() + old_c_size, std::plus<c_data_type>());
			ret.resize(diff + old_c_size);
			//update the indices so that they are the full indices, and not just row indices
			for (auto cur_c_it = ret.index_begin() + old_c_size; cur_c_it < ret.index_end(); ++cur_c_it){
				(*cur_c_it) += this->height()*(b_cur_column + b.width()*cur_tab);
			}
		}
		a_tab_begin = a_temp_end;

	}


    return ret;
//    typedef typename SparseProductReturnType<Derived,otherDerived>::type c_type;
//    typedef typename otherDerived::index_type b_index_type;
//
//    timer_minor=boost::timer::cpu_timer();
//    typename c_type::Indices c_indices;
//    c_indices.reserve(b.dimensionality()*.8); //hard to predict size
//    typename c_type::Data c_data;
//    c_data.reserve(b.dimensionality()*.8); //hard to predict size
//    typename internal::data_type<c_type>::type cur_c_data;
//    std::cout << "Did reserve " << boost::timer::format(timer_minor.elapsed()) << std::endl;
//    timer_minor=boost::timer::cpu_timer();
//    this->print();
//    this->sort(RowMajor);
//    this->print();
//    std::cout << "Did sort " << boost::timer::format(timer_minor.elapsed()) << std::endl;
//    auto a_temp_begin=this->index_begin();
//    auto a_temp_end=a_temp_begin;
//
//    while(a_temp_begin<this->index_end()){
//        auto cur_tab=this->tab(*a_temp_begin);
//        auto cur_row=this->row(*a_temp_begin);
//        a_temp_end=a_temp_begin+1;
//        while(a_temp_end<this->index_end()&&this->tab(*a_temp_end)==cur_tab && this->row(*a_temp_end)==cur_row){
//            a_temp_end++;
//        }
//        for(b_index_type b_columns=0;b_columns<b.width();++b_columns){
//            cur_c_data=0;
//            for(auto a_it=a_temp_begin;a_it<a_temp_end;++a_it){
//                cur_c_data+=this->data_at(a_it)*b(this->column(*a_it),b_columns,cur_tab);
//
//            }
//            if(std::abs(cur_c_data)>SparseTolerance<double>::tolerance){
//                c_data.push_back(cur_c_data);
//                c_indices.push_back(b_columns +(cur_row+cur_tab*this->height())*b.width());
//            }
//        }
//        a_temp_begin=a_temp_end;
//    }
//    c_type ret(std::move(c_data),std::move(c_indices),this->height(),b.width(),this->depth());
//    ret.set_linIdxSequence(RowMajor);
//    std::cout << "Sparse x dense time " << boost::timer::format(timer.elapsed()) << std::endl;
//    this->sort(ColumnMajor);
//    return ret;


}


template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived, otherDerived>::type SparseLatticeBase<Derived>::old_times(SparseLatticeBase<otherDerived> &b){

	this->check_mult_dims(b);
	typedef typename ScalarPromoteType<Derived, otherDerived>::type super_data_type;

	typedef typename SparseProductReturnType<Derived, otherDerived>::type RType;
	typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions
	//iterators to the current tab start and end indexes
	auto a_temp_begin = this->index_begin(), a_temp_end = this->index_begin();
	auto a_index_end = this->index_end();
	auto b_temp_begin = b.index_begin(), b_temp_end = b.index_begin();
	auto b_index_end = b.index_end();
	//maps the rows of each tab of A from compressed indices to full indices
	std::vector<index_type> a_row_map;
	//stores the starting location of each unique column in each tab of A
	std::vector<size_t> a_column_idx;
	//workspace vectors (similar to what's done when CSC matrices are multipled - see Direct Methods for Sparse Linear Systems
	std::vector<size_t> row_marker;
	std::vector<super_data_type> data_collector;
	index_type old_m = this->height();
	index_type cur_tab;
	//initial estimate of size of C
	typename RType::Indices c_indices;

	c_indices.reserve(this->size() + b.size());
	typename RType::Data c_data;
	c_data.reserve(this->size() + b.size());


	this->sort();
	b.sort();

	//determine whether we want to binary search or just scan through elements
	bool a_search_flag = false, b_search_flag = false;
	if (this->depth()*log2(this->size())<this->size())
		a_search_flag = true;
	if (b.depth()*log2(b.size())<b.size())
		b_search_flag = true;

	while (a_temp_begin<a_index_end && b_temp_begin<b_index_end){
		//if we're at the same tab, no work
		if (this->tab(*a_temp_begin) == b.tab(*b_temp_begin)){
			cur_tab = this->tab(*a_temp_begin);
		}
		else if (this->tab(*a_temp_begin)<b.tab(*b_temp_begin)){ //if a's tab is less than b's tab - we need to try to find b's tab in a
			cur_tab = b.tab(*b_temp_begin);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "A-- Tab " << this->tab(*a_temp_begin) << " is less than b tab: " << b.tab(*b_temp_begin) << std::endl;
#endif
			a_temp_begin = this->find_tab_start_idx(cur_tab, a_temp_begin, a_index_end, a_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << " so searched for it and got index " << a_temp_begin - this->index_begin() << std::endl;
#endif
			if (a_temp_begin == a_index_end) //no tab in A is greater than or equal to B's current tab - so we're finished the entire multiplication routine
				break;
			//couldn't find a tab in A equal to B's current tab, but we found one greater than it - so now we need to try to find a matching tab in b
			if (this->tab(*a_temp_begin) != b.tab(*b_temp_begin))
				continue;
		}
		else{
			cur_tab = this->tab(*a_temp_begin);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "B-- Tab " << b.tab(*b_temp_begin) << " is less than A tab: " << this->tab(*a_temp_begin) << std::endl;
#endif
			b_temp_begin = b.find_tab_start_idx(cur_tab, b_temp_begin, b_index_end, b_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << " so searched for it and got index " << b_temp_begin - b.index_begin() << std::endl;
#endif
			if (b_temp_begin == b_index_end) //no tab in B is greater than or equal to A's current tab - so we're finished the entire multiplication routine
				break;
			//couldn't find a tab in B equal to A's current tab, but we found one greater than it - so now we need to try to find a matching tab in B
			if (this->tab(*a_temp_begin) != b.tab(*b_temp_begin))
				continue;
		}

		//find the end of the current tab
		a_temp_end = this->find_tab_end_idx(cur_tab, a_temp_begin, a_index_end, a_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
		std::cout << "A-- Tab " << cur_tab << " begin: " << a_temp_begin - this->index_begin() << " end: " << a_temp_end - this->index_begin() << std::endl;
#endif
		b_temp_end = b.find_tab_end_idx(cur_tab, b_temp_begin, b_index_end, b_search_flag);
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
		std::cout << "B-- Tab " << cur_tab << " begin: " << b_temp_begin - b.index_begin() << " end: " << b_temp_end - b.index_begin() << std::endl;
#endif

		std::sort(a_temp_begin, a_temp_end, [this](index_type lhs, index_type rhs)
		{
			return this->row(lhs)<this->row(rhs);
		});
		//find the number of unique rows
		size_t unique_counter = 1;
		for (auto it = a_temp_begin + 1; it<a_temp_end; ++it)
		{
			if (this->row(*it) != this->row(*(it - 1)))
				unique_counter++;
		}
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
		std::cout << cur_tab << " unique_counter " << unique_counter << std::endl;
#endif
		size_t new_m = unique_counter;
		//resize our row map
		a_row_map.resize(new_m);
		a_row_map.assign(new_m, 0);
		//also resize our workspace
		row_marker.resize(new_m);
		row_marker.assign(new_m, 0);
		data_collector.resize(new_m);
		unique_counter = 0;
		//create a row map from compressed rows to full rows, and also map the current tab's row indices to the compressed form
		index_type old_row = this->row(*a_temp_begin);
		a_row_map[0] = old_row;
		*a_temp_begin = new_m*(this->column(*a_temp_begin) + this->tab(*a_temp_begin)*this->width());
		for (auto it = a_temp_begin + 1; it<a_temp_end; ++it)
		{
			if (old_row != this->row(*it))
			{
				old_row = this->row(*it);
				unique_counter++;
				a_row_map[unique_counter] = old_row;
			}
			*it = unique_counter + new_m*(this->column(*it) + this->tab(*it)*this->width());

		}
		//temporarily set A's height to compressed number of rows
		this->set_height(new_m);
		//sort back to column major
		std::sort(a_temp_begin, a_temp_end);
#ifdef LM_COLUMN_SEARCH
		//find the number of unique columns
		unique_counter = 1;
		for (auto it = a_temp_begin + 1; it<a_temp_end; ++it)
		{
			if (this->column(*it) != this->column(*(it - 1)))
				unique_counter++;
		}

		//resize our column idx marker
		a_column_idx.resize(unique_counter);
		unique_counter = 1;
		//store the location of the start of the first column
		a_column_idx[0] = a_temp_begin - this->index_begin();
		for (auto it = a_temp_begin + 1; it<a_temp_end; ++it)
		{
			if (this->column(*it) != this->column(*(it - 1)))
			{
				//store the where in the index array the new column starts
				a_column_idx[unique_counter++] = it - this->index_begin();
			}
		}
#endif
		//iterate through every element of b
		b_index_type cur_column;
		auto cur_b = b_temp_begin;
		while (cur_b<b_temp_end)
		{
			cur_column = b.column(*cur_b);
			//for each column of b, we start at the beginning of A
			index_iterator a_cur_it = a_temp_begin;
			auto a_cur_column_idx = a_column_idx.begin();
			size_t old_c_size = c_indices.size();
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "Tab " << cur_tab << " Column " << cur_column << ": old_c_size " << old_c_size << std::endl;
#endif
			while (cur_b<b_temp_end && b.column(*cur_b) == cur_column)
			{
				//see Tim Davis' book on Direct Sparse Methods for an explanation of mult_scatter (although this one is slightly different as CSC format isn't used)
				a_cur_it = mult_scatter(a_cur_column_idx, a_column_idx.end(), a_cur_it, a_temp_end, b.row(*cur_b), cur_column, b.data_at(cur_b - b.index_begin()),
					row_marker, data_collector, c_indices);
				cur_b++;
				//if we've reached the end of A's current tab, we're done looking at b's current column
				if (a_cur_it == a_temp_end)
					break;
			}
			//if we've added to c's indices in the current column of B, then clean it up and add to c's data
			//note mult_scatter just pushes the rows of c to c_indices
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "Tab " << cur_tab << " Column " << cur_column << ": new_c_size " << c_indices.size() << std::endl;
#endif
			if (c_indices.size() - old_c_size)
			{
				std::sort(c_indices.begin() + old_c_size, c_indices.end()); //scatter doesn't put the rows in sorted order
				c_data.resize(c_indices.size());
				for (auto it = c_indices.begin() + old_c_size; it<c_indices.end(); ++it)
				{
					c_data[it - c_indices.begin()] = data_collector[(size_t)*it];
					*it = a_row_map[(size_t)*it] + old_m*(cur_column + cur_tab*b.width()); //decompresses the index
				}
			}

		}
		//decompress A's indices
		for (auto it = a_temp_begin; it<a_temp_end; ++it)
		{
			*it = a_row_map[this->row(*it)] + old_m*(this->column(*it) + this->width()*cur_tab);
		}
		//reset A's height to the original height
		this->set_height(old_m);





		a_temp_begin = a_temp_end;
		b_temp_begin = b_temp_end;

	}

	return RType(std::move(c_data), std::move(c_indices), this->height(), b.width(), this->depth());

}




//!Only creates the column index array needed for CSC format
template<class Derived>
void SparseLatticeBase<Derived>::clearCP(){
    this->CP.clear();
	this->CP.shrink_to_fit();
    this->IR.clear();
	this->IR.shrink_to_fit();
}


//!Only creates the column index array needed for CSC format
template<class Derived>
void SparseLatticeBase<Derived>::createCP(index_iterator tab_begin, index_iterator tab_end){

    this->CP.resize(this->width()+1);
    this->IR.resize(tab_end-tab_begin);
    auto cur_column_begin=this->tab(*tab_begin)*this->width()*this->height();
    auto cur_column_end=cur_column_begin+this->height();
    auto cur_it=tab_begin;
    size_t k=0;
    for(index_type i=0;i<this->width();++i){


        this->CP[i] = cur_it-tab_begin;
        while(cur_it<tab_end && *cur_it<cur_column_end){
            this->IR[k++]=*cur_it-cur_column_begin; //we'll use our data array to store just rows
            ++cur_it;
        }

        cur_column_begin=cur_column_end;
        cur_column_end+=this->height();
    }
    this->CP[this->width()] = tab_end-tab_begin;

}

//!Based off of Bulic and Gilbert's DCSC vectors in their hypersparse data structures and their DCSC constructor in CombBlas
template<class Derived>
void SparseLatticeBase<Derived>::createJCCPIR(index_iterator tab_begin, index_iterator tab_end){

        index_type tab_size=tab_end-tab_begin;



        index_type localnzc = 0;
        index_type cur_tab_adder=this->tab(*tab_begin)*this->width()*this->height();        
        auto cur_it=tab_begin;
        //we count columns this way becuase although asymtotically slower, it's faster practically as it avoids divisions to get current column
        
		//need to fix all of these - right now running time is based on width, which is not what we want when we use DCSC
		while(cur_it<tab_end){
			auto cur_column = this->column(*cur_it);
			auto cur_val_end = cur_tab_adder + (cur_column + 1)*this->height();
			++cur_it;
			while (cur_it < tab_end && *cur_it<cur_val_end){
				++cur_it;
            }
			++localnzc;
			
        }

        this->JC.resize(localnzc);
        this->CP.resize(localnzc+1);
        this->IR.resize(this->size());
        index_type jspos = 0;        
        cur_it=tab_begin;        
		
		while (cur_it < tab_end){
			auto cur_column = this->column(*cur_it);
			this->JC[jspos] = cur_column;
			this->CP[jspos++] = cur_it - tab_begin;
			auto cur_val_end = cur_tab_adder + (cur_column+1)*this->height();			
			while (cur_it < tab_end && *cur_it < cur_val_end){
				this->IR[cur_it - tab_begin] = *cur_it - cur_val_end + this->height();
				++cur_it;
			}
		}
		
        this->CP[jspos] = tab_size;


}


//!Based off of Bulic and Gilber's DCSC vectors in their hypersparse data structures and their DCSC constructor in CombBlas. Modified so that chunksize, and hence the divisor, is a power of two (for speed purposes)
template<class Derived>
auto SparseLatticeBase<Derived>::createAUX(index_iterator tab_begin, index_iterator tab_end)->fast_divisor{

        assert(JC.size());
        auto localnzc=JC.size();
        float cf  = static_cast<float>(this->width()+1) / static_cast<float>(localnzc);
        cf=std::floor(log2(cf));
        cf=static_cast<float>(std::pow(2,cf));
        index_type colchunks = static_cast<index_type> ( std::ceil( static_cast<float>(this->width()+1) / std::ceil(cf)) );
        this->AUX.resize(colchunks+1);


        index_type chunksize	= static_cast<index_type>(ceil(cf));
        index_type reg		= 0;
        index_type curchunk	= 0;
        this->AUX[curchunk++] = 0;
        index_type chunk_mult=chunksize;
        for(index_type i = 0; i< (index_type)localnzc; ++i)
        {
            if(this->JC[i] >= chunk_mult)		// beginning of the next chunk
            {
                while(this->JC[i] >= chunk_mult)	// consider any empty chunks
                {
                    this->AUX[curchunk++] = reg;
                    chunk_mult+=chunksize;
                }
            }
            reg = i+1;
        }


        while(curchunk <= colchunks)
        {
            this->AUX[curchunk++] = reg;
        }
        //std::cout << "ChunkSize " << chunksize << std::endl;
        return fast_divisor(chunksize);
}

//!Based off of Bulic and Gilber's DCSC vectors in their hypersparse data structures and their DCSC constructor in CombBlas
template<class Derived>
auto SparseLatticeBase<Derived>::createDCSCVectors(index_iterator tab_begin, index_iterator tab_end)->fast_divisor{

        createJCCPIR(tab_begin,tab_end);
        return createAUX(tab_begin,tab_end);
}


template<class Derived>
void SparseLatticeBase<Derived>::clearDCSCVectors(){

        this->JC.clear();
		this->JC.shrink_to_fit();
        this->CP.clear();
		this->CP.shrink_to_fit();
        this->IR.clear();
		this->IR.shrink_to_fit();
        this->AUX.clear();
		this->AUX.shrink_to_fit();
}


//!should be able to handle all mixtures of hyper-sparse and sparse cases, ultimately this should be altered to use some form of poly-algorithm
template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::operator*(SparseLatticeBase<otherDerived> &b){
    return this->csc_no_accum<true>(b);
}


template <class Derived>
template <bool ColumnSparse,class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::csc_no_accum(SparseLatticeBase<otherDerived> &b){

	using namespace std::chrono;
    this->check_mult_dims(b);
    typedef typename ScalarPromoteType<Derived,otherDerived>::type super_data_type;

    typedef typename SparseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions
    //iterators to the current tab start and end indexes
    



    index_type old_m=this->height();
    index_type cur_tab;
    //initial estimate of size of C
    typename RType::Indices c_indices;

    c_indices.reserve(std::max(this->size(),b.size()));
    typename RType::Data c_data;
    c_data.reserve(c_indices.capacity());


	this->sort(ColumnMajor);
	b.sort(ColumnMajor);

	auto a_temp_begin = this->index_begin(), a_temp_end = this->index_begin();
	auto a_index_end = this->index_end();
	auto b_temp_begin = b.index_begin(), b_temp_end = b.index_begin();
	auto b_index_end = b.index_end();

    //determine whether we want to binary search or just scan through elements
    bool a_search_flag=false, b_search_flag=false;
    if(this->depth()*log2(this->size())<this->size())
        a_search_flag=true;
    if(b.depth()*log2(b.size())<b.size())
        b_search_flag=true;


    while(a_temp_begin<a_index_end && b_temp_begin<b_index_end){


		//if we're at the same tab, no work
        if(this->tab(*a_temp_begin)==b.tab(*b_temp_begin)){
            cur_tab=this->tab(*a_temp_begin);
        }
        else if (this->tab(*a_temp_begin)<b.tab(*b_temp_begin)){ //if a's tab is less than b's tab - we need to try to find b's tab in a
            cur_tab=b.tab(*b_temp_begin);
            a_temp_begin=this->find_tab_start_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
            if(a_temp_begin==a_index_end) //no tab in A is greater than or equal to B's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in A equal to B's current tab, but we found one greater than it - so now we need to try to find a matching tab in b
            if(this->tab(*a_temp_begin)!=b.tab(*b_temp_begin))
                continue;
        }
        else{
            cur_tab=this->tab(*a_temp_begin);
            b_temp_begin=b.find_tab_start_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);
            if(b_temp_begin==b_index_end) //no tab in B is greater than or equal to A's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in B equal to A's current tab, but we found one greater than it - so now we need to try to find a matching tab in B
            if(this->tab(*a_temp_begin)!=b.tab(*b_temp_begin))
                continue;
        }


        a_temp_end=this->find_tab_end_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
        b_temp_end=b.find_tab_end_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);

        MultHelper<Derived,ColumnSparse> multHelper(*this);
        multHelper.initMultHelpers(a_temp_begin,a_temp_end);
		

        //iterate through every element of b
        b_index_type cur_column;

        auto cur_tab_adder=cur_tab*this->height()*b.width();
        b_index_type cur_b_tab_adder=cur_tab*b.height()*b.width();
        auto cur_b = b_temp_begin;

        while (cur_b<b_temp_end)
        {

            //for each column of b, we start at the beginning of A
            size_t old_c_size=c_indices.size();

            cur_column = b.column(*cur_b);
			auto cur_offset_end=cur_b_tab_adder+(cur_column+1)*b.height();
			//for each column of b, we start at the beginning of A

#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "Tab " << cur_tab << " Column " << cur_column << ": old_c_size " << old_c_size << " *cur_b " << *cur_b << " cur_offset_end " << cur_offset_end << std::endl;
#endif
            while (cur_b<b_temp_end && *cur_b<cur_offset_end){
                //see Tim Davis' book on Direct Sparse Methods for an explanation of mult_scatter (although this one is slightly different as the accumulator isn't used)
                mult_scatter(a_temp_begin,b.row(*cur_b),cur_column,b.data_at(cur_b),c_indices,c_data,multHelper);
                cur_b++;
            }

            if(c_indices.size()-old_c_size )
            {

                if(c_indices.size()-old_c_size>1){
                    internal::RadixSortInPlace_PowerOf2Radix_Unsigned(c_indices.begin()+old_c_size, c_data.begin()+old_c_size,c_indices.size()-old_c_size,this->height() );

                    auto diff=collect_duplicates_function(c_indices.begin()+old_c_size, c_indices.end(), c_data.begin()+old_c_size,std::plus<super_data_type>());
                    c_indices.resize(diff+old_c_size);
                    c_data.resize(diff+old_c_size);
                }

                cur_offset_end=cur_tab_adder+this->height()*cur_column;
                for(auto it=c_indices.begin()+old_c_size; it<c_indices.end(); ++it)                {

					*it +=cur_offset_end; //decompresses the index
                }
            }
        }
		
        //high_resolution_clock::time_point t2 = high_resolution_clock::now();
        //std::cout << "Pure mult time: " << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;
        //std::cout << "Pure Mult Time:" << boost::timer::format(hash_t.elapsed()) << std::endl;
    }

    return RType(std::move(c_data),std::move(c_indices),this->height(),b.width(),this->depth(),true);

}


template<class Derived>
template<class b_index_type, class b_data_type, class ret_index_type, class ret_data_type,bool RowSparse>
void SparseLatticeBase<Derived>::mult_scatter(const index_iterator a_begin, const b_index_type b_row, const b_index_type b_column, const b_data_type beta, std::vector<ret_index_type> & c_indices,std::vector<ret_data_type> & c_data,const MultHelper<Derived,RowSparse> &multHelper) const
{





    index_type a_column_begin,a_column_end;
    multHelper.getColumnOffset(b_row,a_column_begin,a_column_end);

    if(a_column_begin==a_column_end)
        return;

	auto _begin=this->IR.begin()+a_column_begin;
	auto _end=this->IR.begin()+a_column_end;
	auto _data_it=this->data_begin()+(a_begin-this->index_begin())+a_column_begin;

    if(c_indices.size()+(_end-_begin)>c_indices.capacity()){
        c_indices.reserve(2*c_indices.capacity());
        c_data.reserve(2*c_data.capacity());
    }

	c_data.resize(c_data.size() + (_end - _begin));
	for (int i = c_data.size() - (_end - _begin); i < c_data.size(); ++i){
		c_data[i] = beta*(*_data_it++);
	}
	c_indices.resize(c_indices.size() + (_end - _begin));
	std::copy(_begin, _end, c_indices.end() - (_end - _begin));

	//for (auto cur_row = _begin; cur_row<_end; ++cur_row){
	//	c_indices.push_back(*cur_row); //push back the compressed row (this will have to be remapped afterwards)
	//	c_data.push_back(beta*(*_data_it++));
	//}


}


template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::outer_times(SparseLatticeBase<otherDerived> &b){


    this->check_mult_dims(b);
    typedef typename ScalarPromoteType<Derived,otherDerived>::type super_data_type;

    typedef typename SparseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions
    //iterators to the current tab start and end indexes
    auto a_temp_begin=this->index_begin(), a_temp_end=this->index_begin();
    auto a_index_end=this->index_end();
    auto b_temp_begin=b.index_begin(),b_temp_end=b.index_begin();
    auto b_index_end=b.index_end();



    index_type old_m=this->height();
    index_type cur_tab;
    //initial estimate of size of C
    typename RType::Indices c_indices;

    c_indices.reserve(10000);
    typename RType::Data c_data;
    c_data.reserve(c_indices.capacity());


	this->sort(ColumnMajor);
	b.sort(RowMajor);




    //determine whether we want to binary search or just scan through elements
    bool a_search_flag=false, b_search_flag=false;
    if(this->depth()*log2(this->size())<this->size())
        a_search_flag=true;
    if(b.depth()*log2(b.size())<b.size())
        b_search_flag=true;


    while(a_temp_begin<a_index_end && b_temp_begin<b_index_end){


		//if we're at the same tab, no work
        if(this->tab(*a_temp_begin)==b.tab(*b_temp_begin)){
            cur_tab=this->tab(*a_temp_begin);
        }
        else if (this->tab(*a_temp_begin)<b.tab(*b_temp_begin)){ //if a's tab is less than b's tab - we need to try to find b's tab in a
            cur_tab=b.tab(*b_temp_begin);
            a_temp_begin=this->find_tab_start_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
            if(a_temp_begin==a_index_end) //no tab in A is greater than or equal to B's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in A equal to B's current tab, but we found one greater than it - so now we need to try to find a matching tab in b
            if(this->tab(*a_temp_begin)!=b.tab(*b_temp_begin))
                continue;
        }
        else{
            cur_tab=this->tab(*a_temp_begin);
            b_temp_begin=b.find_tab_start_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);
            if(b_temp_begin==b_index_end) //no tab in B is greater than or equal to A's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in B equal to A's current tab, but we found one greater than it - so now we need to try to find a matching tab in B
            if(this->tab(*a_temp_begin)!=b.tab(*b_temp_begin))
                continue;
        }
        a_temp_end=this->find_tab_end_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
        b_temp_end=b.find_tab_end_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);



        auto cur_tab_adder=cur_tab*this->height()*b.width();
        auto tab_size=this->height()*b.width();
        b_index_type b_cur_row;
        index_type a_cur_column;
        auto a_cur_it=a_temp_begin;
        auto b_cur_it=b_temp_begin;
        auto old_c_size=c_indices.size();
        while(a_cur_it<a_temp_end && b_cur_it <b_temp_end){
            a_cur_column=this->column(*a_cur_it);
            b_cur_row=b.row(*b_cur_it);
            //std::cout << "a_cur_column " << a_cur_column << " b_cur_row " << b_cur_row << std::endl;
            if(a_cur_column<b_cur_row){
				++a_cur_it;
				while (a_cur_it < a_temp_end && this->column(*a_cur_it) < b_cur_row)
					++a_cur_it;
            }
            else if(b_cur_row<a_cur_column){
				++b_cur_it;
				while (b_cur_it < b_temp_end && b.row(*b_cur_it) < a_cur_column)
					++b_cur_it;
            }
            else{
                auto b_cur_end=b_cur_it;
                auto a_cur_end=a_cur_it;

                while(a_cur_end<a_temp_end && this->column(*a_cur_end)==a_cur_column){
                    a_cur_end++;
                    //std::cout << "a looking: " << this->row(*a_cur_end) << " " << this->column(*(a_cur_end++)) << std::endl;
                }

                while(b_cur_end<b_temp_end && b.row(*b_cur_end)==a_cur_column){
                    //std::cout << "b looking: " << b.row(*b_cur_end) << " " << b.column(*(b_cur_end++)) << std::endl;
                    b_cur_end++;
                }
               // std::cout << " a size " << (a_cur_end-a_cur_it) << " b size " << (b_cur_end-b_cur_it) << std::endl;
                if(c_indices.capacity()<c_indices.size()+(a_cur_end-a_cur_it)*(b_cur_end-b_cur_it)){
                    c_indices.reserve(2*c_indices.capacity());
                    c_data.reserve(2*c_data.capacity());
                }
                for(;b_cur_it<b_cur_end;++b_cur_it){



                    for(auto a_cur_inner=a_cur_it;a_cur_inner<a_cur_end;++a_cur_inner){
                        c_indices.push_back(this->row(*a_cur_inner)+b.column(*b_cur_it)*old_m);
                        c_data.push_back(this->data_at(a_cur_inner)*b.data_at(b_cur_it));
                    }

                }
                a_cur_it=a_cur_end;
            }



        }
		//std::cout << "Outer f " << c_indices.size() - old_c_size << std::endl;
        if(c_indices.size()-old_c_size){
            internal::RadixSortInPlace_PowerOf2Radix_Unsigned(c_indices.begin()+old_c_size, c_data.begin()+old_c_size,c_indices.size()-old_c_size,tab_size );
    //internal::Introsort(c_indices.begin()+old_c_size,c_indices.end(),std::less<index_type>(),internal::DualSwapper<index_iterator,data_iterator>(c_indices.begin()+old_c_size,c_data.begin()+old_c_size));
            auto diff=collect_duplicates_function(c_indices.begin()+old_c_size, c_indices.end(), c_data.begin()+old_c_size,std::plus<super_data_type>());
            c_indices.resize(diff+old_c_size);
            c_data.resize(diff+old_c_size);
            for(auto it=c_indices.begin()+old_c_size;it<c_indices.end();++it)
                *it+=cur_tab_adder;

        }

        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;
    }

        //std::cout << "bucket count " << _bucket_count << " load factor " << _load_factor << " size " << _size << std::endl;







    return RType(std::move(c_data),std::move(c_indices),this->height(),b.width(),this->depth(),true);

}


template<class Derived>
template<class b_index_type, class b_data_type,class ret_index_type,class super_data_type,bool ColSparse>
bool SparseLatticeBase<Derived>::mult_scatter(index_type a_tab_offset,b_index_type b_row,b_index_type b_column,b_data_type beta,
	std::vector<b_index_type>& row_marker, std::vector<super_data_type> & data_collector, std::vector<ret_index_type> & c_indices, const MultHelper<Derived, ColSparse> & multHelper)
{

    index_type a_column_begin,a_column_end;
    multHelper.getColumnOffset(b_row,a_column_begin,a_column_end);

	if (a_column_begin == a_column_end)
		return false;
    //std::cout << "Index of column at " << b_row << " is " << a_cur_it-this->index_begin() << std::endl;

    if(c_indices.capacity()<c_indices.size()+a_column_end-a_column_begin)
        c_indices.reserve(c_indices.capacity()*2);

	auto _begin=this->IR.begin()+a_column_begin;
	auto _end=this->IR.begin()+a_column_end;
	auto _data_it = this->data_begin() + a_tab_offset + a_column_begin;

	bool need_sorted = false;
    for(auto cur_row=_begin;cur_row<_end;++cur_row){
        if ((b_index_type)(row_marker[*cur_row])<b_column+1){
            row_marker[*cur_row]=b_column+1;
            c_indices.push_back(*cur_row); //push back the compressed row (this will have to be remapped afterwards)
            //std::cout << "Beta: " << beta << " a_data " << this->data_at(a_cur_it-this->index_begin()) << " index " << a_cur_it-this->index_begin() << " result " << beta*this->data_at(a_cur_it-this->index_begin()) << std::endl;
            data_collector[*cur_row]=beta*(*_data_it++);
			need_sorted = true;
        }
        else{
            //std::cout << "Beta: " << beta << " a_data " << this->data_at(a_cur_it-this->index_begin()) << " index " << a_cur_it-this->index_begin() << " result " << beta*this->data_at(a_cur_it-this->index_begin()) << std::endl;
            data_collector[*cur_row]+=beta*(*_data_it++);
        }
    }
	return need_sorted;


}


//#define LM_SPARSE_LATTICE_MULT_DEBUG
template <class Derived>
template <bool ColSparse,class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::csc_times(SparseLatticeBase<otherDerived> &b)
{


    this->check_mult_dims(b);

    typedef typename ScalarPromoteType<Derived,otherDerived>::type super_data_type;
    typedef typename SparseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions


    sort(ColumnMajor); //tab/column major for A
    b.sort(ColumnMajor); //tab/column major for B
    std::vector<b_index_type> row_marker;
    std::vector<super_data_type> data_collector(this->height());
    //iterators for indices and data

    auto a_temp_begin=this->index_begin();
    auto a_index_end=this->index_end();

    auto b_temp_begin=b.index_begin();
    auto b_index_end=b.index_end();

    decltype(b_temp_begin) b_temp_end;
    decltype(a_temp_begin) a_temp_end;

    typename RType::Indices c_indices;
    typename RType::Data c_data;
	c_indices.reserve(std::max(this->size(), b.size()));
    c_data.reserve(c_indices.capacity());


    //determine whether we want to binary search or just scan through elements
    bool a_search_flag=false, b_search_flag=false;
    if(this->depth()*log2(this->size())<this->size())
        a_search_flag=true;
    if(b.depth()*log2(b.size())<b.size())
        b_search_flag=true;
    index_type cur_tab;
    while(a_temp_begin<a_index_end && b_temp_begin<b_index_end){
        row_marker.assign(this->height(),0);
        //if we're at the same tab, no work
        if(this->tab(*a_temp_begin)==b.tab(*b_temp_begin)){
            cur_tab=this->tab(*a_temp_begin);
        }
        else if (this->tab(*a_temp_begin)<b.tab(*b_temp_begin)){ //if a's tab is less than b's tab - we need to try to find b's tab in a
            cur_tab=b.tab(*b_temp_begin);

            a_temp_begin=this->find_tab_start_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);

            if(a_temp_begin==a_index_end) //no tab in A is greater than or equal to B's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in A equal to B's current tab, but we found one greater than it - so now we need to try to find a matching tab in b
            if(this->tab(*a_temp_begin)!=cur_tab)
                continue;
        }
        else{
            cur_tab=this->tab(*a_temp_begin);

            b_temp_begin=b.find_tab_start_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);

            if(b_temp_begin==b_index_end) //no tab in B is greater than or equal to A's current tab - so we're finished the entire multiplication routine
                break;
            //couldn't find a tab in B equal to A's current tab, but we found one greater than it - so now we need to try to find a matching tab in B
            if(b.tab(*b_temp_begin)!=cur_tab)
                continue;
        }

        //find the end of the current tab
        a_temp_end=this->find_tab_end_idx(cur_tab,a_temp_begin,a_index_end,a_search_flag);
        b_temp_end=b.find_tab_end_idx(cur_tab,b_temp_begin,b_index_end,b_search_flag);

        MultHelper<Derived,ColSparse> multHelper(*this);
        multHelper.initMultHelpers(a_temp_begin,a_temp_end);

        b_index_type cur_column;
        b_index_type cur_tab_b_adder=cur_tab*b.height()*b.width();
        auto cur_tab_adder=cur_tab*this->height()*b.width();
		auto cur_b = b_temp_begin;
		auto cur_a_tab_offset = a_temp_begin - this->index_begin();
		while (cur_b<b_temp_end)
		{
			cur_column = b.column(*cur_b);
			auto cur_offset_end=cur_tab_b_adder+(cur_column+1)*b.height();
			
			size_t old_c_size = c_indices.size();
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "Tab " << cur_tab << " Column " << cur_column << ": old_c_size " << old_c_size << " *cur_b " << *cur_b << " cur_offset_end " << cur_offset_end << std::endl;
#endif
			auto beta = b.data_at(cur_b - b.index_begin());
			auto test = b.row(*cur_b);

			//do the first entry in the current b column, as it can be done faster
			index_type a_column_begin, a_column_end;			
			multHelper.getColumnOffset(b.row(*cur_b), a_column_begin, a_column_end);
			if (a_column_begin == a_column_end){
				cur_b++;
				continue;
			}

			if (c_indices.capacity() < c_indices.size() + a_column_end - a_column_begin){
				c_indices.reserve(c_indices.capacity() * 2);
			}
			c_indices.resize(c_indices.size() + a_column_end - a_column_begin);
			auto _begin = this->IR.begin() + a_column_begin;
			auto _end = this->IR.begin() + a_column_end;
			std::copy(_begin, _end, c_indices.begin() + old_c_size);

			
			auto _data_it = this->data_begin() + cur_a_tab_offset + a_column_begin;
			
			for (auto cur_row = _begin; cur_row<_end; ++cur_row){				
				row_marker[*cur_row] = cur_column + 1;						
				data_collector[*cur_row] = beta*(*_data_it++);				
			}	
			
				
			
			cur_b++;

			//now do the rest of the rows
			bool need_sorted = false;
			while (cur_b<b_temp_end && *cur_b<cur_offset_end)
			{
				beta = b.data_at(cur_b - b.index_begin());
				//see Tim Davis' book on Direct Sparse Methods for an explanation of mult_scatter (although this one is slightly different as CSC format isn't used)
				need_sorted = need_sorted | mult_scatter(cur_a_tab_offset, b.row(*cur_b), cur_column, beta, row_marker, data_collector, c_indices, multHelper);
				cur_b++;
			}
			//if we've added to c's indices in the current column of B, then clean it up and add to c's data
			//note mult_scatter just pushes the rows of c to c_indices
#ifdef LM_SPARSE_LATTICE_MULT_DEBUG
			std::cout << "Tab " << cur_tab << " Column " << cur_column << ": new_c_size " << c_indices.size() << std::endl;
#endif
			
			if (need_sorted){
				std::sort(c_indices.begin() + old_c_size, c_indices.end()); //scatter doesn't put the rows in sorted order
			}
				
			cur_offset_end = cur_tab_adder + this->height()*cur_column; //indices currently only store the rows, so convert them to full indices based on current column and current tab

			
			c_data.reserve(c_indices.capacity());
			c_data.resize(c_indices.size());
			for (auto it = c_indices.begin() + old_c_size; it < c_indices.end(); ++it)
			{
				c_data[it - c_indices.begin()] = data_collector[*it];
				*it += cur_offset_end;

			}
				
			

		}

        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;


    }
    //std::cout << "************Finished MULT************" << std::endl;
    return RType(std::move(c_data),std::move(c_indices),this->height(),b.width(),this->depth());
}

//template <class Derived>
//inline void SparseLatticeBase<Derived>::to_matrix(const std::vector<index_type> &row_map, const std::vector<index_type> &inner_indices,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_cm& mat) const
//{
//
//
//
//
//
//    index_type outer_val=0;
//    auto mat_rows=mat.innerIndexPtr();
//    auto mat_columns=mat.outerIndexPtr();
//    typename std::vector<index_type>::const_iterator cur_row;
//    for(auto cur_it=inner_indices.begin(); cur_it<inner_indices.end(); cur_it++)
//    {
//
//        //std::cout << "Cur inner index " << *cur_it << "column of lattice " << column(index(*t_begin)) << std::endl;
//        *mat_columns=outer_val;
//        //std::cout << "column set to " << *mat_columns << std::endl;
//        mat_columns++;
//        cur_row=row_map.begin();
//        while (t_end-t_begin>0 && column(index_val(*t_begin))<*cur_it){
//
//            t_begin++;
//        }
//
//        //less than greater than operators don't seem to work with tupleit (id3eally this function should just use index/data iterators instead.
//        while (t_end-t_begin>0 && column(index_val(*t_begin))==*cur_it)
//        {
//            cur_row=std::lower_bound(cur_row,row_map.end(),row(index_val(*t_begin)));
//            *mat_rows=cur_row-row_map.begin();
//            //std::cout << "seting row to " << *mat_rows << std::endl;
//            mat_rows++;
//            outer_val++;
//            t_begin++;
//
//
//        }
//        if(t_end-t_begin<1){
//
//            break;
//        }
//
//
//
//    }
//
//    while(mat_columns<=mat.outerIndexPtr()+mat.outerSize()){
//        //std::cout << "final column " << mat_columns-mat.outerIndexPtr() << "set to " << mat.nonZeros() << std::endl;
//        *mat_columns=mat.nonZeros();
//        mat_columns++;
//    }
//
//
//}
//
//
//template <class Derived>
//inline void SparseLatticeBase<Derived>::to_matrix_rowmajor(const std::vector<index_type> &inner_indices,const std::vector<index_type> &column_map,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_rm& mat) const
//{
//
//
//
//
//    index_type outer_val=0;
//    auto mat_columns=mat.innerIndexPtr();
//    auto mat_rows=mat.outerIndexPtr();
//    typename std::vector<index_type>::const_iterator cur_column;
//    for(auto cur_it=inner_indices.begin(); cur_it<inner_indices.end(); cur_it++)
//    {
//
//        //std::cout << "Cur inner index " << *cur_it << "row of lattice " << row(index(*t_begin)) << std::endl;
//        *mat_rows=outer_val;
//        //std::cout << "row set to " << *mat_rows << std::endl;
//        mat_rows++;
//        cur_column=column_map.begin();
//
//
//        while (t_end-t_begin>0 && row(index_val(*t_begin))<*cur_it)
//            t_begin++;
//
//        while (t_end-t_begin>0 && row(index_val(*t_begin))==*cur_it)
//        {
//            cur_column=std::lower_bound(cur_column,column_map.end(),column(index_val(*t_begin)));
//            *mat_columns=cur_column-column_map.begin();
//            //std::cout << "seting column to " << *mat_columns << std::endl;
//            outer_val++;
//            mat_columns++;
//            t_begin++;
//        }
//        if(t_end-t_begin<1)
//            break;
//
//
//    }
//    while(mat_rows<=mat.outerIndexPtr()+mat.outerSize()){
//        //std::cout << "final row " << mat_rows-mat.outerIndexPtr() << "set to " << mat.nonZeros() << std::endl;
//        *mat_rows=mat.nonZeros();
//        mat_rows++;
//    }
//
//
//}


namespace{
//!Helper class to pull least squares or householder inversion
template<bool LSQR,class EigenSparseMatrix>
struct sparse_lattice_solver;

template<class EigenSparseMatrix>
struct sparse_lattice_solver<true,EigenSparseMatrix> {

    typedef Eigen::SparseQR<EigenSparseMatrix,Eigen::COLAMDOrdering<typename EigenSparseMatrix::Index>> solverType;

    solverType _solver;
    void initialize_solver(const EigenSparseMatrix & _matrix,SolveInfo &_solveInfo){




        _solver.compute(_matrix);
        if(_solver.info()!=Eigen::Success){
            _solveInfo=RankDeficient;

        }
        else{
            _solveInfo=FullyRanked;
        }
    }

    template<typename BType, typename CType>
    void _solve(const BType & b, CType & c)
    {
        _solver._solve(b,c);
    }
    auto info()->decltype(_solver.info())
    {
        return _solver.info();
    }

    auto lastErrorMessage()->decltype(_solver.lastErrorMessage()){
        return _solver.lastErrorMessage();
    }

//        std::cout << "analyze succeeded " << std::endl;
//        _solver.factorize(_matrix);
//        if(_solver.info()!=Eigen::Success){
//            std::cout << "factoring failed " << std::endl;
//            _solveInfo=RankDeficient;
//        }




};

template<class EigenSparseMatrix>
struct sparse_lattice_solver<false,EigenSparseMatrix> {

    typedef Eigen::SparseLU<EigenSparseMatrix> solverTypeLU;

    //typedef Eigen::SimplicialLDLT<EigenSparseMatrix> solverTypeChol;
    solverTypeLU _solverLU;
    //solverTypeChol _solverChol;

    void initialize_solver(const EigenSparseMatrix & _matrix,SolveInfo &_solveInfo){


        _solverLU.analyzePattern(_matrix);
        _solverLU.factorize(_matrix);

		if (_solverLU.info() != Eigen::Success){

            _solveInfo=RankDeficient;
        }
        else{

            _solveInfo=LeastSquares;
        }


    }

    template<typename BType, typename CType>
    void _solve(const BType & b, CType & c){

            _solverLU._solve(b,c);


    }
    auto info()->decltype(_solverLU.info()){

            return _solverLU.info();
    }

    auto lastErrorMessage()->decltype(_solverLU.lastErrorMessage()){

            return _solverLU.lastErrorMessage();


    }

};


//template<class EigenSparseMatrix>
//struct sparse_lattice_solver<false, EigenSparseMatrix> {
//
//	typedef Eigen::SimplicialLDLT<EigenSparseMatrix> solverTypeLU;
//
//	//typedef Eigen::SimplicialLDLT<EigenSparseMatrix> solverTypeChol;
//	solverTypeLU _solverLU;
//	//solverTypeChol _solverChol;
//
//	void initialize_solver(const EigenSparseMatrix & _matrix, SolveInfo &_solveInfo){
//
//
//		_solverLU.compute(_matrix);
//		if (_solverLU.info() != Eigen::Success){
//			_solveInfo = RankDeficient;
//
//		}
//		else{
//			_solveInfo = FullyRanked;
//		}
//
//
//	}
//
//	template<typename BType, typename CType>
//	void _solve(const BType & b, CType & c){
//
//		_solverLU._solve(b, c);
//
//
//	}
//	auto info()->decltype(_solverLU.info()){
//
//		return _solverLU.info();
//	}
//
//	std::string lastErrorMessage(){
//
//		return "";
//
//
//	}
//
//};
}

template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(SparseLatticeBase<otherDerived> &b)
{


    this->check_solve_dims(b);






    if (this->width()==this->height()){
        return perform_solve<otherDerived,false>(b);
    }
    else if(this->width()<this->height()){
        return perform_solve<otherDerived,true>(b);

    }
    else
        throw LatticeParameterException("Only square or overdetermined systems are supported");







}

template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::lsqr_solve(SparseLatticeBase<otherDerived> &b)
{


    this->check_solve_dims(b);






    if (this->width()<this->height()){
        return perform_solve<otherDerived,true>(b);
    }
//    else if(this->width()<this->height()){
//        return perform_solve<otherDerived,true>(b);
//
//    }
    else
        throw LatticeParameterException("Only over-determined systems are supported");







}

template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(const DenseLatticeBase<otherDerived> &b)
{


    this->check_solve_dims(b);






    if (this->width()==this->height()){
        return perform_solve<otherDerived,false>(b);
    }
    else if (this->width()<this->height()){
        return perform_solve<otherDerived,true>(b);

    }
    else
        throw LatticeParameterException("Only square or overdetermined systems are supported");






}

template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::lsqr_solve(const DenseLatticeBase<otherDerived> &b)
{


    this->check_solve_dims(b);






    if (this->width()<this->height()){
        return perform_solve<otherDerived,true>(b);
    }

    else
        throw LatticeParameterException("Only over-determined systems are supported");






}

template <class Derived>
template <class otherDerived,bool LSQR>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::perform_solve(SparseLatticeBase<otherDerived> &b)
{


    //std::cout << "Entered sparse solve sparse " << std::endl;
    this->check_solve_dims(b);



    typedef typename SparseSolveReturnType<Derived,otherDerived>::type c_type;
    typedef typename c_type::vector_type c_vector_type;




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

    std::vector<typename MappedSparseMatrix_cm::Index> a_rows; //we make it MappedSparseMatrix_cm::Index, incase index_type differs from MappedSparseMatrix_cm::Index






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
			//SparseMatrix_cm A = A_temp;
            //get solver
			typedef sparse_lattice_solver<LSQR, MappedSparseMatrix_cm> solverType;
            solverType _solver;
            _solver.initialize_solver(A,this->mSolveInfo);

            if(this->mSolveInfo==RankDeficient)
            {
                std::stringstream t;
                t << "Could not perform facotrization on tab " << k << " due to error : " << _solver.lastErrorMessage() << ".";
                //std::cout << A << std::endl;

                throw RankDeficientException(t.str());
            }

            //loop through every column of b
            while(b_temp_begin<b_temp_end){
                b_vector.setZero(); //reset dense temp vector

                auto b_cur_column=b.column(*b_temp_begin);
                while(b_temp_begin<b_temp_end && b.column(*b_temp_begin)==b_cur_column){
                    b_vector(b.row(*b_temp_begin))=this->convert(b.derived().data_at(b_temp_begin-b.index_begin()));
                    b_temp_begin++;
                }


                //get wrapper for corresponding column of lattice c
                c_vector_type c_vector=c.column_vector(b_cur_column,k);
                //solve and store in lattice c
                _solver._solve(b_vector,c_vector);
                if(_solver.info()!=Eigen::Success)
                {
                    this->mSolveInfo=RankDeficient;
                    std::stringstream t;
                    t << "Solution process on tab " << k << " and column "<< b_cur_column << "of RHS failed due to error : " << _solver.lastErrorMessage() << ".";
                    throw RankDeficientException(t.str());
                }


            }

        }
        else
        {
            this->mSolveInfo=RankDeficient;
            std::stringstream t;
            t << "Rank deficient tab. Tab " << k << "has zero entries.";
            throw RankDeficientException(t.str());

        }
        //update current tab location
        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;

    }


    //std::cout << "Exited sparse solve sparse " << std::endl;
    return c;
}
template <class Derived>
template <class otherDerived,bool LSQR>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::perform_solve(const DenseLatticeBase<otherDerived> &b)
{

#ifdef LM_SPARSE_LATTICE_SOLVE_DEBUG
	std::cout << "Entered sparse solve with dense " << std::endl;
#endif
	

    this->check_solve_dims(b);
    typedef typename SparseSolveReturnType<Derived,otherDerived>::type c_type;





    sort(ColumnMajor); //tab/column major for A

#ifdef LM_SPARSE_LATTICE_SOLVE_DEBUG
    std::cout << "Finished sort" << std::endl;
#endif
    //iterators for for indices and data
    index_iterator a_temp_begin=this->index_begin();
    auto a_index_end=this->index_end();
    index_iterator a_temp_end;


    //hold compressed indices for Sparse matrix of each tab
    std::vector<index_type> a_columns;
    //since all columns must have at least one non-zero entry to be fully-ranked, we can just do direct mapping to CSC
    a_columns.resize(this->width()+1);

    typedef const typename  Eigen::MappedSparseMatrix<data_type,Eigen::ColMajor,index_type> MappedSparseMatrix_cm; //Sparse Matrix type for A


    c_type c(this->width(),b.width(),this->depth());   //create dense lattice to return and allocate memory
#ifdef LM_SPARSE_LATTICE_SOLVE_DEBUG
    std::cout << "Allocated" << std::endl;
#endif
    for (int k=0; k<this->depth(); k++) //loop through every tab
    {


        a_temp_begin=this->find_tab_start_idx(k,a_temp_begin,a_index_end);

        //find the last occurence of current tab in A lattices
        a_temp_end=this->find_tab_end_idx(k,a_temp_begin,a_index_end);

        //tab from a must have nonzeros
        if (a_temp_end!=a_temp_begin)
        {


            //now we temporarily remap a_indices to rows
            *a_temp_begin=this->row(*a_temp_begin); //set first to 0 row
            auto a_cur_it=a_temp_begin+1;
            size_t cur_column=0; //if size_t can't hold the number of columns, we're in big trouble anyway
            a_columns[0]=0;
            while(a_cur_it<a_temp_end){
                if(this->column(*a_cur_it)!=cur_column){
                    //check to make sure we didn't skip a column - if we did then we have a rank-deficient tab b/c of a zero-column
                    if(this->column(*a_cur_it)-cur_column>1){
                        std::stringstream t;
                        t << "Rank deficient tab. Column " << cur_column << " in tab " << k << " of LHS has zero entries.";
                        throw RankDeficientException(t.str());
                    }
                    cur_column++;
#ifdef LM_SPARSE_LATTICE_SOLVE_DEBUG
                    std::cout << "Moving to column " << cur_column << std::endl;
#endif
                    a_columns[cur_column]=a_cur_it-a_temp_begin;
                }
                *a_cur_it=this->row(*a_cur_it); //remap indices to rows
                ++a_cur_it;
            }

            //end value of outer index array
            a_columns.back()=a_temp_end-a_temp_begin;

            //created a CCS matrix by mapping row and column vectors and also the pre-existing data of *this lattice
            MappedSparseMatrix_cm A=MappedSparseMatrix_cm(this->height(),this->width(),a_columns.back(),&a_columns[0],&(*a_temp_begin),&(*(data_begin()+(a_temp_begin-this->index_begin())))); //map data to a compressed column matrix

#ifdef LM_SPARSE_LATTICE_SOLVE_DEBUG
			std::cout << "Made ccs matrix" << std::endl;
#endif
            //compute LU decomposition
            //get solver
            typedef sparse_lattice_solver<LSQR,MappedSparseMatrix_cm> solverType;
            solverType _solver;
            _solver.initialize_solver(A,this->mSolveInfo);
            if(this->mSolveInfo==RankDeficient)
            {
                std::stringstream t;
                t << "Could not perform facotrization on tab " << k << " due to error : " << _solver.lastErrorMessage() << ".";
                throw RankDeficientException(t.str());
            }


            //for some reason, the QR solver throws an exception if it's passed two EigenMap matrices. But if doesn't when passed with individual EigenMap columns
            for(size_t col_idx=0;col_idx<(size_t)b.width();++col_idx){
                auto c_vector=c.column_vector(col_idx,k);
                auto b_vector=b.column_vector(col_idx,k);
                _solver._solve(b_vector,c_vector);
               
                if(_solver.info()!=Eigen::Success)
                {
                    this->mSolveInfo=RankDeficient;
                    std::stringstream t;
                    t << "Solution process on tab " << k << "  col " << col_idx << " failed due to error : " << _solver.lastErrorMessage() << ".";
                    throw RankDeficientException(t.str());
                }

            }



             //now map back to A's full index, using the temp CSC matrix that was created
             //***ACtually, this should happen even if the solution process failed - must change in all sparse routines!
            for(auto column_it=a_columns.begin();column_it<a_columns.end()-1;++column_it){
                auto cur_column=column_it-a_columns.begin();
                for(auto row_it=a_temp_begin+*column_it;row_it<a_temp_begin+*(column_it+1);++row_it)
                    *row_it=*row_it+this->height()*(cur_column+this->width()*k);
            }


        }
        else
        {
            this->mSolveInfo=RankDeficient;
            std::stringstream t;
            t << "Rank deficient tab. Tab " << k << "has zero entries.";
            throw RankDeficientException(t.str());

        }
        //update current tab location
        a_temp_begin=a_temp_end;


    }

#ifdef LM_SPARSE_LATTICE_SOLVE_DEBUG
    std::cout << "Finished sparse solve with dense " << std::endl;
#endif
    return c;
}


//template <class Derived>
//template <class otherDerived>
//typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(SparseLatticeBase<otherDerived> &b){
//    //use delayed parsing, so the static assert will only trigger if the function is actually used within a compilation unit
//    struct fake : std::false_type{};
//    static_assert(fake::value,"You must define LIBMIA_USE_SPARSE_SOLVE (or turn it on if using CMake) and build and link to SuperLU to perform sparse solution of equations.");
//}
//template <class Derived>
//template <class otherDerived>
//typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(const DenseLatticeBase<otherDerived> &b){
//    //use delayed parsing, so the static assert will only trigger if the function is actually used within a compilation unit
//    struct fake : std::false_type{};
//    static_assert(fake::value,"You must define LIBMIA_USE_SPARSE_SOLVE (or turn it on if using CMake) and build and link to SuperLU to perform sparse solution of equations.");
//}





/*! @} */




} //libMIA





#endif // SPARSELATTICEBASE_H
