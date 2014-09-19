// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef DENSEMIA_H_INCLUDED
#define DENSEMIA_H_INCLUDED

#include <type_traits>
#include <iostream>
#include <algorithm>


#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>


#include "LibMIAException.h"
#include "LibMIAUtil.h"
#include "IndexUtil.h"
#include "DenseMIABase.h"
#include "ImplicitMIA.h"


//\defgroup
namespace LibMIA
{

/** \addtogroup mia Multi-Index Array Classes
 *  @{
 */
namespace internal
{

template<typename T,size_t _order>
struct data_type<DenseMIA<T,_order> >
{
    typedef T type;
};

template<typename T,size_t _order>
struct data_type_ref<DenseMIA<T,_order> >
{
    typedef T & type;
};
template<typename T,size_t _order>
struct const_data_type_ref<DenseMIA<T,_order> >
{
    typedef const T & type;
};

template<typename T,size_t _order>
struct index_type<DenseMIA<T,_order> >
{
    typedef long long type;
};

template<typename T,size_t _order>
struct order<DenseMIA<T,_order> >
{
    constexpr static size_t value=_order;
};




template<typename T,size_t _order>
struct data_iterator<DenseMIA<T,_order> >
{
    typedef T* restrict_libmia type;
};

template<typename T,size_t _order>
struct const_data_iterator<DenseMIA<T,_order> >
{
    typedef const T* restrict_libmia type;
};

template<typename T,size_t _order>
struct FinalDerived<DenseMIA<T,_order> >
{
    typedef DenseMIA<T,_order> type;
};

}




//!  MIA class for dense data.
/*!
  Supports addition, multiplication, and solution of, possibly over-determined, systems of
  linear equations. By default MIA will own underlying raw data, which can be reallocated and resized.

  \tparam T   the datatype of individual elements.
  \tparam _order   the order (number of indices) of the MIA.
*/
template <class T, size_t _order>
class DenseMIA: public DenseMIABase<DenseMIA<T,_order> >
{





public:

    //! raw data_type
    typedef typename internal::data_type<DenseMIA>::type data_type;
    //! raw data_type ref
    typedef typename internal::data_type_ref<DenseMIA>::type data_type_ref;
    //! raw data_type ref
    typedef typename internal::const_data_type_ref<DenseMIA>::type const_data_type_ref;
    //! raw index_type
    typedef typename internal::index_type<DenseMIA>::type index_type;


    //! iterator type for iterating directly through raw data
    typedef typename internal::data_iterator<DenseMIA>::type data_iterator;
    //! iterator type for iterating directly through raw data
    typedef typename internal::const_data_iterator<DenseMIA>::type const_data_iterator;
    //! raw data pointer type
    typedef T* raw_pointer;
    //! order of the MIA
    constexpr static size_t mOrder=_order;
    //! smart pointer type used to reference raw data
    typedef std::unique_ptr<T restrict_libmia [] > smart_raw_pointer;


    //! final derived type
    typedef typename internal::FinalDerived<DenseMIA>::type FinalDerived;

    typedef typename DenseMIABase<DenseMIA>::accumulator_type accumulator_type;
    typedef typename DenseMIABase<DenseMIA>::fast_accumulator_type fast_accumulator_type;
    typedef typename DenseMIABase<DenseMIA>::multiplier_type multiplier_type;

private:

    smart_raw_pointer m_smart_raw_ptr;

    bool hasOwnership;


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

    struct null_deleter
    {
        void operator()(void const *) const
        {
        }
    };

    //!  Constructs empty DenseMIA
    DenseMIA():DenseMIABase<DenseMIA<T,_order> >(),m_smart_raw_ptr(nullptr),hasOwnership(true)
    {
    }

    //!  Constructs DenseMIA of specified size with a given raw data pointer.
    /*!
        \tparam[in] _dims The dimensions size of data. Will assert a compile failure is size is different than _order
        \tparam[in] scalar_data Raw pointer of scalar data
        \tparam[in] _ownership If false, DenseMIA will not own data and caller is responsible for deleting scalar_data. Also DenseMIA will not be
                                allowed to resize scalar_data, but can assign individual entries. If true, scalar data must be allocated using new [], otherwise
                                behaviour is undefined upon destruction of DenseMIA object
    */
    template<class array_index_type>
    DenseMIA(const std::array<array_index_type,_order> &_dims,T* scalar_data,bool _ownership=true):DenseMIABase<DenseMIA<T,_order> >(_dims),hasOwnership(_ownership)
    {

        m_smart_raw_ptr.reset(scalar_data);



    }




    template<class array_index_type>//, class otherDerived>
    DenseMIA(const std::array<array_index_type,_order> &_dims,DenseLattice<data_type>&& _lat,bool _ownership=true):DenseMIA(_dims,_lat.release_memptr(),_ownership)
    {
        this->mSolveInfo=_lat.solveInfo();
    }

    //!  Copy constructor.
    /*!


    */
    DenseMIA(const DenseMIA& restrict_libmia otherMIA)  :DenseMIABase<DenseMIA<T,_order> >(otherMIA.dims()),hasOwnership(true)
    {

        this->mSolveInfo=otherMIA.solveInfo();
        m_smart_raw_ptr.reset(new T[this->dimensionality()]);


        if(this->dimensionality()>=PARALLEL_TOL){
            #pragma omp parallel for
            for(size_t idx=0;idx<this->dimensionality();++idx)
            {
                this->atIdx(idx)=otherMIA.atIdx(idx);
            }
        }
        else
        {
            for(size_t idx=0;idx<this->dimensionality();++idx)
            {
                this->atIdx(idx)=otherMIA.atIdx(idx);
            }
        }


    }

    //!  Move constructor.
    /*!


    */
    DenseMIA(DenseMIA&& otherMIA):DenseMIABase<DenseMIA<T,_order> >(otherMIA.dims()),m_smart_raw_ptr(nullptr),hasOwnership(otherMIA.ownsData())
    {




        this->mSolveInfo=otherMIA.solveInfo();
        m_smart_raw_ptr.swap(otherMIA.m_smart_raw_ptr);



    }


    //!  Copy constructor.
    /*!
        If otherMIA's datatype is different than this->data_type, then individual entries will be converted.

    */
    template<class otherMIAType, typename boost::enable_if<typename internal::is_DenseMIA<otherMIAType>::type,int>::type=0>
    DenseMIA(const otherMIAType& restrict_libmia otherMIA):DenseMIABase<DenseMIA<T,_order> >(otherMIA.dims()),hasOwnership(true)
    {


        this->mSolveInfo=otherMIA.solveInfo();
        static_assert(internal::order<otherMIAType>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");
        m_smart_raw_ptr.reset(new T[this->dimensionality()]);


        if(this->dimensionality()>=PARALLEL_TOL){
        #pragma omp parallel for
            for(size_t idx=0;idx<this->dimensionality();++idx)
            {
                this->atIdx(idx)=this->convert(otherMIA.atIdx(idx));
            }
        }
        else
        {
            for(size_t idx=0;idx<this->dimensionality();++idx)
            {
                this->atIdx(idx)=this->convert(otherMIA.atIdx(idx));
            }
        }




    }

    //!  Copy constructor for SparseMIA.
    /*!
        If otherMIA's datatype is different than this->data_type, then individual entries will be converted.

    */
    template<class otherMIAType, typename boost::enable_if<typename internal::is_SparseMIA<otherMIAType>::type,int>::type=0>
    DenseMIA(const otherMIAType& restrict_libmia otherMIA):DenseMIABase<DenseMIA<T,_order> >(otherMIA.dims()),hasOwnership(true)
    {



        static_assert(internal::order<otherMIAType>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");
        this->straight_assign(otherMIA);






    }

    //!clears the raw data and sets dimensions to zero
    void clear(){
        if (ownsData()){
            auto _newDims=this->dims();
            std::fill(_newDims.begin(),_newDims.end(),0);
            this->set_dims(_newDims);

            m_smart_raw_ptr.reset(nullptr);

        }
        else
            throw MIAMemoryException("Cannot clear data if DenseMIA does not own data.");


    }


    //!  Constructs DenseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<typename... Dims>
    DenseMIA(Dims... dims):DenseMIABase<DenseMIA<T,_order> > {dims...}, m_smart_raw_ptr(new T[this->dimensionality()]),hasOwnership(true)
    {



        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        this->zeros();


    }

    //!  Constructs DenseMIA of specified size.
    /*!
        Scalar data will be set to zero

        \param[in] dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than _order

    */
    template<class array_index_type>
    DenseMIA(const std::array<array_index_type,mOrder>& _dims):DenseMIABase<DenseMIA<T,_order> > (_dims), m_smart_raw_ptr(new T[this->dimensionality()]),hasOwnership(true)
    {


        this->zeros();

    }

    //!  Constructs DenseMIA of specified size with a given raw data pointer.
    /*!

        \param[in] scalar_data Raw pointer of scalar data
        \param[in] _ownership If false, DenseMIA will not own data and caller is responsible for deleting scalar_data. Also DenseMIA will not be
                                allowed to resize scalar_data, but can assign individual entries. If true, scalar data must be allocated using new [], otherwise
                                behaviour is undefined upon destruction of DenseMIA object
        \param[in] _dims variadic parameter to specify size. Will assert a compile failure is number of parameters is different than mOrder or if datatype of variables making
                            up _dims are not convertible to index_type
    */
    template<typename... Dims>
    DenseMIA(T* scalar_data,bool _ownership,Dims... dims):DenseMIABase<DenseMIA<T,_order> > {dims...},hasOwnership(_ownership)
    {

        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");

        m_smart_raw_ptr.reset(scalar_data);



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


    //!  Assignment based on given order.
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up

        \param[in] otherMIA the other MIA
        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {2,0,1} this->at(x,y,z)==otherMIA.at(y,z,x).

    */
    template<typename otherDerived,typename index_param_type>
    void assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order);

    template<typename other_data_type,typename index_param_type>
    void assign(const ImplicitMIA<other_data_type,mOrder>& restrict_libmia otherMIA,const std::array<index_param_type,_order>& index_order);

    //!  Assignment based on given order using an rvalue reference.
    /*!

        If the data_type of otherMIA is not the same as this, the scalar data will be converted. The function allows a user to specify
        a permutation of indices to shuffle around the scalar data. Will assert compile failure if the orders of the two MIAs don't match up

        \param[in] otherMIA the other MIA
        \param[in] index_order The assignment order, given for otherMIA. E.g., if order is {2,0,1} this->at(x,y,z)==otherMIA.at(y,z,x).

        If *this owns its data, an in-place permutation of otherMIA's data is performed, and then *this is set to refer to otherMIA's raw data

    */
    template<typename index_param_type>
    void assign(DenseMIA&& otherMIA,const std::array<index_param_type,_order>& index_order);

    //!  Assignment operator.
    /*!
        \param[in] otherMIA
        If the data_type of otherMIA is not the same as this, the scalar data will be converted.
    */
    template<typename other_data_type>
    DenseMIA& operator=(const DenseMIA<other_data_type,mOrder>& otherMIA);

    template<typename other_data_type>
    DenseMIA& operator=(const ImplicitMIA<other_data_type,mOrder>& restrict_libmia otherMIA);

    template<typename otherDerived>
    DenseMIA& operator=(const SparseMIABase<otherDerived>& otherMIA);

    DenseMIA& operator=(const DenseMIA& otherMIA);

//    //!  Straight-out assignment.
//    DenseMIA& operator=(const DenseMIA& otherMIA){
//        if(this!=&otherMIA)
//            return *this=static_cast<const DenseMIABase<DenseMIA>&>(otherMIA);
//        else
//            return *this;
//    }
    //!Move assignment
    DenseMIA& operator=(DenseMIA&& restrict_libmia otherMIA) restrict_libmia
    {




        if(this==&otherMIA)
            return *this;

        if (this->ownsData()){ //we can only swap if data is owned by *this


            //std::cout << "Better be here " << std::endl;
            m_smart_raw_ptr.swap(otherMIA.m_smart_raw_ptr);



            std::swap(this->hasOwnership,otherMIA.hasOwnership);
            this->set_dims(otherMIA.dims());
            this->mSolveInfo=otherMIA.solveInfo();

            return *this;
        }
        else{ //perform copy assignment

            return *this=static_cast<DenseMIA&>(otherMIA);
        }

    }



    bool ownsData() const{
        return hasOwnership;
    }

    //! Returns a raw pointer to the scalar data
    inline T* restrict_libmia raw_data_ptr() const{
        return m_smart_raw_ptr.get();
    }

    //! Returns a raw pointer to the scalar data and releases ownership - caller must deallocate data using delete[]
    inline T* restrict_libmia release_raw_data() {
        hasOwnership=false;
        return raw_data_ptr();
    }

    //! Iterator to the beginning of the raw data
    inline data_iterator data_begin()
    {
        return raw_data_ptr();

    }

    //! Iterator to the end of the raw data
    inline data_iterator data_end()
    {
        return raw_data_ptr()+size();

    }

        //! Returns scalar data at given linear index
    inline const_data_type_ref atIdx(index_type idx) const{

        //return lin index
        return *(data_begin()+idx);
    }

    //! Returns scalar data at given linear index
    inline data_type_ref atIdx(index_type idx){

        //return lin index
        return *(data_begin()+idx);
    }

    //! Iterator to the beginning of the raw data
    inline const_data_iterator data_begin() const
    {
        return raw_data_ptr();

    }

    //! Iterator to the end of the raw data
    inline const_data_iterator data_end() const
    {
        return raw_data_ptr()+size();

    }

    //! Returns size of raw data. For dense cases, this is the same as dimensionality
    std::size_t size() const
    {

        return this->dimensionality();

    }



    //! If this has ownership, the raw data will be deallocated
    ~DenseMIA(){
        if(!hasOwnership)
            m_smart_raw_ptr.release();
    }




    //! Performs destructive add (+=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is added like: this->at(x,y,z)+=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,typename index_param_type>
    DenseMIA & plus_equal(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){

        std::plus<data_type> op;
        merge(b.derived(),op,index_order);
        return *this;
    }


    //! Performs destructive subtract (-=).
    /*!
        \param[in] index_order The assignment order, given for b. E.g., if order is {3,1,2}, then each data_value is subtracted like: this->at(x,y,z)-=b.at(z,x,y).
        Will assert a compile failure is size of index_order is not the same as this->mOrder
    */
    template<class otherDerived,typename index_param_type>
    DenseMIA & minus_equal(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order){
        std::minus<data_type> op;
        merge(b.derived(),op,index_order);
        return *this;
    }

    template<typename index_param_type>
    void inplace_permute(const std::array<index_param_type,internal::order<DenseMIA>::value> & reshuffle_order);




    //!
    /*!
        Based on An Optimal Index Reshuffle Algorithm for Multidimensional Arrays and Its Applications for Parallel Architectures by Chris H.Q. Ding and the modification
        in Sec. III A by Jie et al.'s article: A High Efficient In-place Transposition Scheme for Multidimensional Arrays
    */
    template<typename... Dims>
    void inplace_permute(Dims... dims){
        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");

        std::array<index_type, mOrder> reshuffle_order{{dims...}};
        inplace_permute(reshuffle_order);

    }




//    //!
//    /*!
//        An Idea tried to speed-up inplace_permutation, but ended up being slower
//    */
//    template<typename... Dims>
//    void inplace_permute_new(Dims... dims){
//        static_assert(internal::check_mia_constructor<DenseMIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
//
//        std::array<index_type, mOrder> reshuffle_order{{dims...}};
//        inplace_permute_new(reshuffle_order);
//
//    }
//  void inplace_permute_new(const std::array<index_type,mOrder> & reshuffle_order);



    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
    template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    void    merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,_order>& index_order);
    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
    template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_SparseMIA<otherDerived>, int >::type = 0>
    void  merge(const MIA<otherDerived> &b,const Op& op,const std::array<index_param_type,_order>& index_order);

    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
    template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type = 0>
    void    merge(const MIA<otherDerived> &b,const Op& op);
    //! Common routine for merge operations, such as add or subtract. Templated on the Op binary operator.
    template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_SparseMIA<otherDerived>, int >::type = 0>
    void  merge(const MIA<otherDerived> &b,const Op& op);

    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom.
    /*!
        For DenseMIAs, this function calls toLatticeCopy
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R_size, size_t C_size, size_t T_size>
    DenseLattice<data_type> toLatticeExpression(const std::array<idx_typeR,R_size> & row_indices, const std::array<idx_typeC,C_size> & column_indices,const std::array<idx_typeT,T_size> & tab_indices) const
    {
        return this->toLatticeCopy(row_indices, column_indices, tab_indices);

    }

    //! Flattens the MIA to a Lattice by permuting the data in-place.
    /*!
        \param[in] row_indices indices to map to the lattice rows - will perserve ordering
        \param[in] column_indices indices to map to the lattice columns - will perserve ordering
        \param[in] tab_indices indices to map to the lattice tabs - will perserve ordering
        \return MappedDenseLattice class that wraps this's raw data
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R_size, size_t C_size, size_t T_size>
    MappedDenseLattice<data_type> toLatticePermute(const std::array<idx_typeR,R_size> & row_indices, const std::array<idx_typeC,C_size> & column_indices,const std::array<idx_typeT,T_size> & tab_indices);


    //! Flattens the MIA to a Lattice. This function is called in MIA expressions by MIA_Atom when the MIA in question is a temp object.
    /*!
        For DenseMIAs, this function calls toLatticePermute
    */
    template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R_size, size_t C_size, size_t T_size>
    auto toLatticeDiscard(const std::array<idx_typeR,R_size> & row_indices, const std::array<idx_typeC,C_size> & column_indices,const std::array<idx_typeT,T_size> & tab_indices)
    ->decltype(this->toLatticePermute(row_indices, column_indices, tab_indices))
    {
        return this->toLatticePermute(row_indices, column_indices, tab_indices);

    }

    auto toStraightLattice(size_t number_of_row_indices, size_t number_of_column_indices) const ->MappedDenseLattice<data_type>;

protected:
    template<typename other_data_type>
    DenseMIA& straight_assign(const DenseMIA<other_data_type,mOrder>& otherMIA);

    template<typename otherDerived>
    DenseMIA& straight_assign(const SparseMIABase<otherDerived>& otherMIA);



private:


    template <class D1,class D2,bool D3,size_t D4> friend class MIA_Atom;
    template <class E1> friend class DenseMIABase;
    template <class F1> friend class MIA;



};



template<class T, size_t _order>
template<typename index_param_type>
void DenseMIA<T,_order>::inplace_permute(const std::array<index_param_type,internal::order<DenseMIA>::value>& reshuffle_order){

    //first check that the reshuffle_order array isn't just {0,1,...mOrder}
    //if it is, we do nothing
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

    size_t check;
    for(check=0;check<reshuffle_order.size();++check){
        if(reshuffle_order[check]!=(index_param_type)check)
            break;
    }
    if (check==reshuffle_order.size())
        return;

    boost::dynamic_bitset<> bit_array(this->dimensionality()); // all 0's by default
    std::array<index_type,mOrder> new_dims; //stores new dimensions
    internal::reorder_from(this->dims(),reshuffle_order,new_dims); //get new dims
    index_type ioffset;

    const auto dim_accumulator=internal::createDimAccumulator(new_dims,reshuffle_order);
    const auto multiplier=internal::createMultiplier(this->dims());

    //create a function that converts from a linIdx of the permuted array to a linIdx of the old array (see Jie et al.'s article: A High Efficient In-place Transposition Scheme for Multidimensional Arrays)
    //EDIT - using the lambda slowed down the implementation, so it's just hard-coded now
//    auto func=[this,&dim_accumulator](const index_type from_lin_idx){
//        index_type to_lin_idx=0;
//        index_type multiplier=1;
//        for(size_t i=0;i<mOrder;++i){
//            to_lin_idx+=(from_lin_idx/dim_accumulator[i])%this->dim(i)*multiplier; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly
//            multiplier*=this->dim(i);
//        }
//        return to_lin_idx;
//    };
    index_type touch_ctr=0;
    //iterate through the entire bit array
    for(index_type start_idx=0;start_idx<this->dimensionality();++start_idx){
        //if we've found a location that hasn't been touched, we've found a new vacancy cycle
        if(bit_array[start_idx]==0){

            auto temp=this->atIdx(start_idx); //get the start of the cycle
            ioffset=start_idx; //permuted array linIdx
            while(true){

                size_t ioffset_next=0;

                for(size_t i=0;i<mOrder;++i){
                    ioffset_next+=((size_t)ioffset/dim_accumulator[i])%(size_t)this->dim(i)*multiplier[i]; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly

                }
                //ioffset_next=func(ioffset); //get the location of ioffset in the old array


                bit_array[ioffset]=1; //touch the current data location
                ++touch_ctr;
                if(ioffset_next==start_idx){ //if we've cycled to the start, then finish the cycle and break
                    if(ioffset_next!=ioffset)
                        this->atIdx(ioffset)=temp;
                    break;
                }

                this->atIdx(ioffset)=this->atIdx(ioffset_next); //set the data element in the permuted array to its location in the old array
                ioffset=ioffset_next; //update the location in the current cycle

            }
            if (touch_ctr==this->dimensionality()){
                break;
            }

        }


    }
    this->set_dims(new_dims);


}



//template<typename Derived>
//void DenseMIABase<Derived>::inplace_permute_new(const std::array<index_type,mOrder>& reshuffle_order){
//    //boost::dynamic_bitset<> bit_array(this->dimensionality()); // all 0's by default
//    std::array<index_type,mOrder> new_dims; //stores new dimensions
//    internal::reorder_from(this->dims(),reshuffle_order,new_dims); //get new dims
//    index_type ioffset_next,ioffset;
//    std::string dummy;
//    std::array<index_type,mOrder> dim_accumulator; //precompute the demoninators needed to convert from linIdx to a full index, using new_dims
//    for(size_t i=0;i<mOrder;++i){
//        dim_accumulator[i]=std::accumulate(new_dims.begin(),new_dims.begin()+i,1,std::multiplies<index_type>());
//    }
//
//    dim_accumulator=internal::reorder_to(dim_accumulator,reshuffle_order); //reorder the denominators based on the reshuffle order
//    //create a function that converts from a linIdx of the permuted array to a linIdx of the old array (see Jie et al.'s article: A High Efficient In-place Transposition Scheme for Multidimensional Arrays)
//    auto func=[this,&dim_accumulator](const index_type from_lin_idx){
//        index_type to_lin_idx=0;
//        index_type multiplier=1;
//        for(size_t i=0;i<mOrder;++i){
//            to_lin_idx+=(from_lin_idx/dim_accumulator[i])%this->dim(i)*multiplier; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly
//            multiplier*=this->dim(i);
//        }
//        return to_lin_idx;
//    };
//    index_type touch_ctr=0;
//    //iterate through the entire bit array
//    index_type bottom_idx=-1;
//    index_type start_idx;
//    bool do_bottom=true;
//    index_type top_idx=this->dimensionality();
//    while(true){
//    //for(index_type start_idx=0;start_idx<this->dimensionality();++start_idx){
//        //if we've found a location that hasn't been touched, we've found a new vacancy cycle
//
//
//        if (do_bottom){
//            bottom_idx++;
//            start_idx=bottom_idx;
//
//        }
//        else{
//            top_idx--;
//            start_idx=top_idx;
//        }
//        //std::cout <<"Start idx " << start_idx << std::endl;
//        ioffset_next=func(start_idx); //get the location of ioffset in the old array
//        //std::cout <<"ioffset_next " << ioffset_next << std::endl;
//        while(ioffset_next>bottom_idx&& ioffset_next<top_idx){
//            index_type old_one=ioffset_next;
//            ioffset_next=0;
//            index_type multiplier=1;
//            for(size_t i=0;i<mOrder;++i){
//                ioffset_next+=(old_one/dim_accumulator[i])%this->dim(i)*multiplier; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly
//                multiplier*=this->dim(i);
//            }
//            //ioffset_next=func(ioffset_next);
//            //std::cout <<"ioffset_next " << ioffset_next << std::endl;
//        }
//        //std::cin >> dummy;
//
//        do_bottom=!do_bottom;
//        //std::cout <<"Past" << std::endl;
//        if(ioffset_next==start_idx){ //found a new cycle
//            //std::cout <<"In" << std::endl;
//            auto temp=this->atIdx(start_idx); //get the start of the cycle
//            ioffset=start_idx; //permuted array linIdx
//            while(true){
//                //std::cout << "ioffset " << ioffset << std::endl;
//                index_type ioffset_next=0;
//                index_type multiplier=1;
//                for(size_t i=0;i<mOrder;++i){
//                    ioffset_next+=(ioffset/dim_accumulator[i])%this->dim(i)*multiplier; //use the shuffled denominators to compute shuffle full indices, then convert linIdx on the fly
//                    multiplier*=this->dim(i);
//                }
//
//                //ioffset_next=func(ioffset); //get the location of ioffset in the old array
//                //std::cout << "ioffset next" << ioffset_next << std::endl;
//                //std::cin >> dummy;
//                ++touch_ctr;
//                if(ioffset_next==start_idx){ //if we've cycled to the start, then finish the cycle and break
//                    if(ioffset_next!=ioffset)
//                        this->atIdx(ioffset)=temp;
//                    break;
//                }
//
//                this->atIdx(ioffset)=this->atIdx(ioffset_next); //set the data element in the permuted array to its location in the old array
//                ioffset=ioffset_next; //update the location in the current cycle
//            }
//            if (touch_ctr==this->dimensionality()){
//                break;
//            }
//
//        }
//
//
//    }
//    this->m_dims=new_dims;
//
//}


template<class T, size_t _order>
template<typename other_data_type>
DenseMIA<T,_order>& DenseMIA<T,_order>::operator=(const ImplicitMIA<other_data_type,mOrder>& otherMIA){
    //if the otherMIA is implicit, we need to get its explicit values in case it's part of an expression that includes *this


    if(!hasOwnership && this->m_dimensionality!=otherMIA.dimensionality()){
        throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    auto temp=otherMIA.template make_explicit<data_type>();
    return *this=std::move(temp); //if *this owns its data then no copy will take place, just a rewrap of raw memory



}

template<class T, size_t _order>
template<typename other_data_type>
DenseMIA<T,_order>& DenseMIA<T,_order>::operator=(const DenseMIA<other_data_type,mOrder>& otherMIA)
{

    return this->straight_assign(otherMIA);

}
template<class T, size_t _order>
DenseMIA<T,_order>& DenseMIA<T,_order>::operator=(const DenseMIA<T,_order>& otherMIA)
{

    return this->straight_assign(otherMIA);

}


template<class T, size_t _order>
template<typename otherDerived>
DenseMIA<T,_order>& DenseMIA<T,_order>::operator=(const SparseMIABase<otherDerived>& otherMIA)
{

    static_assert(internal::order<otherDerived>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");
    return this->straight_assign(otherMIA);

}


template<class T, size_t _order>
template<typename other_data_type>
DenseMIA<T,_order>& DenseMIA<T,_order>::straight_assign(const DenseMIA<other_data_type,mOrder>& otherMIA)
{




    if(this->dimensionality()!=otherMIA.dimensionality()){
        if(hasOwnership){
            this->set_dims(otherMIA.dims());
            smart_raw_pointer temp_ptr(new T[this->dimensionality()]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();

        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    else if(this->dims()!=otherMIA.dims())
        this->set_dims(otherMIA.dims());

    this->mSolveInfo=otherMIA.solveInfo();

    //std::cout << "about to convert" << std::endl;
    for(size_t idx=0;idx<this->dimensionality();++idx){
        //std::cout << idx << std::endl;
        this->atIdx(idx)=this->convert(otherMIA.atIdx(idx));
    }

    return *this;

}


template<class T, size_t _order>
template<typename otherDerived>
DenseMIA<T,_order>& DenseMIA<T,_order>::straight_assign(const SparseMIABase<otherDerived>& otherMIA){

    static_assert(internal::order<otherDerived>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");

    if(this->dimensionality()!=otherMIA.dimensionality()){
        if(hasOwnership){

            this->set_dims(otherMIA.dims());
            smart_raw_pointer temp_ptr(new T[this->dimensionality()]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();

        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    else if(this->dims()!=otherMIA.dims()){
        this->set_dims(otherMIA.dims());
    }
    this->zeros();
    this->mSolveInfo=otherMIA.solveInfo();

    //std::cout << "about to convert" << std::endl;
    for(auto it=otherMIA.index_begin();it<otherMIA.index_end();++it){

        //std::cout << idx << std::endl;
        this->atIdx(otherMIA.convert_to_default_linIdxSequence(*it))=this->convert(otherMIA.data_at(it));
    }

    return *this;

}



template<class T, size_t _order>
template<typename other_data_type,typename index_param_type>
void DenseMIA<T,_order>::assign(const ImplicitMIA<other_data_type,mOrder>& otherMIA,const std::array<index_param_type,_order>& index_order)
{


    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

    //if the otherMIA is implicit, we need to get its explicit values in it refers to an expression that *this (eg a=b+a)
    if(!hasOwnership && this->m_dimensionality!=otherMIA.dimensionality()){
        throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    *this=otherMIA.template make_explicit<data_type>(index_order);


//    print_array(otherMIA.dims(),"otherMIA.dims()");
//    print_array(index_order,"index_order");
//    print_array(this->m_dims,"this->m_dims");




}



template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void DenseMIA<T,_order>::assign(const SparseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{


    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::order<otherDerived>::value==mOrder,"Order of MIAs must be the same to perform copy construction.");

    if(this->dimensionality()!=otherMIA.dimensionality()){
        if(hasOwnership){
            auto new_dims=this->dims();
            internal::reorder_from(otherMIA.dims(),index_order,new_dims);
            this->set_dims(new_dims);

            smart_raw_pointer temp_ptr(new T[this->dimensionality()]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();

        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
     else{
        auto new_dims=otherMIA.dims();
        internal::reorder_from(otherMIA.dims(),index_order,new_dims);
        this->set_dims(new_dims);

    }

    this->zeros();
    this->mSolveInfo=otherMIA.solveInfo();
    index_type curIdx=0;
    auto lhsOrder=internal::reverseOrder(index_order);
    accumulator_type dim_accumulator;
    fast_accumulator_type fast_dim_accumulator;
    multiplier_type multiplier;
    internal::create_shuffle_needs(otherMIA.dims(),this->dims(),lhsOrder,dim_accumulator,fast_dim_accumulator,multiplier);


    //print_array(dim_accumulator,"dim_accumulator");
    //print_array(index_order,"index_order");



    for(auto it=otherMIA.index_begin();it<otherMIA.index_end();++it){

        //std::cout << idx << std::endl;
        auto default_idx=otherMIA.convert_to_default_linIdxSequence(*it);
        auto shuffled_idx=internal::reShuffleLinearIndex(default_idx,multiplier,fast_dim_accumulator,dim_accumulator);
        this->atIdx(shuffled_idx)=this->convert(otherMIA.data_at(it));
    }


}


template<class T, size_t _order>
template<typename otherDerived,typename index_param_type>
void DenseMIA<T,_order>::assign(const DenseMIABase<otherDerived>& otherMIA,const std::array<index_param_type,_order>& index_order)
{


    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

    if(boost::is_same<DenseMIA, otherDerived>::value && this==&otherMIA){
        //std::cout << "NO! " << std::endl;
        this->inplace_permute(index_order);
        return;
    }

    if(this->dimensionality()!=otherMIA.dimensionality()){
        if(hasOwnership){
            auto new_dims=otherMIA.dims();
            internal::reorder_from(otherMIA.dims(),index_order,new_dims);
            this->set_dims(new_dims);
            smart_raw_pointer temp_ptr(new T[this->dimensionality()]);
            m_smart_raw_ptr.swap(temp_ptr);
            temp_ptr.release();

        }
        else
            throw new MIAMemoryException("Cannot assign to MIA that doesn't own underlying data if dimensionality is different");
    }
    else{
        auto new_dims=otherMIA.dims();
        internal::reorder_from(otherMIA.dims(),index_order,new_dims);
        this->set_dims(new_dims);
    }

    this->mSolveInfo=otherMIA.solveInfo();
    index_type curIdx=0;


    accumulator_type dim_accumulator;
    fast_accumulator_type fast_dim_accumulator;
    multiplier_type multiplier;
    internal::create_shuffle_needs(this->dims(),otherMIA.dims(),index_order,dim_accumulator,fast_dim_accumulator,multiplier);

    for(auto this_it=this->data_begin(); this_it<this->data_end(); ++this_it)
    {

        *this_it=this->convert(otherMIA.atIdx(internal::reShuffleLinearIndex(curIdx++,multiplier,fast_dim_accumulator,dim_accumulator)));

    }


}

template<class T, size_t _order>
template<typename index_param_type>
void DenseMIA<T,_order>::assign(DenseMIA<T,_order>&& otherMIA,const std::array<index_param_type,_order>& index_order)
{

    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");

    if(this==&otherMIA){
        this->inplace_permute(index_order);
        return;
    }
    if(!this->ownsData()){
        this->assign(static_cast<DenseMIA&>(otherMIA), index_order); //copy if *this doesn't own its data
        return;
    }
    //otherwise we are free to swap, so do it, after a permute of course.
    otherMIA.inplace_permute(index_order);
    m_smart_raw_ptr.swap(otherMIA.m_smart_raw_ptr);

    this->hasOwnership=otherMIA.ownsData();
    this->set_dims(otherMIA.dims());
    this->mSolveInfo=otherMIA.solveInfo();



}

template<class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
void  DenseMIA<T,_order>::merge(const MIA<otherDerived>  & restrict_libmia b,const Op& op,const std::array<index_param_type,_order>& index_order) restrict_libmia
{

    #ifdef LIBMIA_CHECK_DIMS
    this->check_merge_dims(b,index_order);
    #endif
    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");



    accumulator_type dim_accumulator;
    fast_accumulator_type fast_dim_accumulator;
    multiplier_type multiplier;
    internal::create_shuffle_needs(this->dims(),b.dims(),index_order,dim_accumulator,fast_dim_accumulator,multiplier);

    if(this->dimensionality()>=PARALLEL_TOL){
        #pragma omp parallel for
        for(index_type curIdx=0; curIdx<this->dimensionality(); ++curIdx)
        {
            this->atIdx(curIdx)=op(this->atIdx(curIdx),this->convert(b.atIdx(internal::reShuffleLinearIndex(curIdx,multiplier,fast_dim_accumulator,dim_accumulator))));


        }
    }
    else{
       for(index_type curIdx=0; curIdx<this->dimensionality(); ++curIdx)
        {
            this->atIdx(curIdx)=op(this->atIdx(curIdx),this->convert(b.atIdx(internal::reShuffleLinearIndex(curIdx,multiplier,fast_dim_accumulator,dim_accumulator))));

        }


    }



}

template<class T, size_t _order>
template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_DenseMIA<otherDerived>, int >::type>
void  DenseMIA<T,_order>::merge(const MIA<otherDerived>  & restrict_libmia b,const Op& op) restrict_libmia
{

    #ifdef LIBMIA_CHECK_DIMS
    this->check_merge_dims(b);
    #endif







    if(this->dimensionality()>=PARALLEL_TOL){
        #pragma omp parallel for
        for(index_type curIdx=0; curIdx<this->dimensionality(); ++curIdx)
        {
            this->atIdx(curIdx)=op(this->atIdx(curIdx),this->convert(b.atIdx(curIdx)));


        }
    }
    else{
       for(index_type curIdx=0; curIdx<this->dimensionality(); ++curIdx)
        {
            this->atIdx(curIdx)=op(this->atIdx(curIdx),this->convert(b.atIdx(curIdx)));

        }


    }



}


template<class T, size_t _order>
auto DenseMIA<T,_order>::toStraightLattice(size_t number_of_row_indices, size_t number_of_column_indices) const ->MappedDenseLattice<data_type>
{




    index_type row_size=std::accumulate(this->dims().begin(),this->dims().begin()+number_of_row_indices,1,std::multiplies<index_type>());
    index_type column_size=std::accumulate(this->dims().begin()+number_of_row_indices,this->dims().begin()+number_of_row_indices+number_of_column_indices,1,std::multiplies<index_type>());
    index_type tab_size=std::accumulate(this->dims().begin()+number_of_row_indices+number_of_column_indices,this->dims().end(),1,std::multiplies<index_type>());

    return MappedDenseLattice<data_type>(raw_data_ptr(),row_size, column_size, tab_size);






}

template<class T, size_t _order>
template< class idx_typeR, class idx_typeC, class idx_typeT, size_t R_size, size_t C_size, size_t T_size>
auto DenseMIA<T,_order>::toLatticePermute(const std::array<idx_typeR,R_size> & row_indices, const std::array<idx_typeC,C_size> & column_indices,const std::array<idx_typeT,T_size> & tab_indices) ->MappedDenseLattice<data_type>
{

    static_assert(internal::check_index_compatibility<index_type,idx_typeR>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeC>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeT>::type::value,"Must use an array convertable to index_type");
    static_assert(R_size+C_size+T_size==mOrder,"Size of all three arrays must equal mOrder");
    //statically check number of indices match up
    size_t row_size=1, column_size=1, tab_size=1;
    std::array<index_type,R_size+C_size+T_size> permute_order;
    internal::concat_arrays(row_indices,column_indices,tab_indices,permute_order);

    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    this->inplace_permute(permute_order);
    for(size_t i=0;i<R_size;++i)
        row_size*=this->dim(i);
    for(size_t i=R_size;i<C_size+R_size;++i)
        column_size*=this->dim(i);
    for(size_t i=C_size+R_size;i<mOrder;++i)
        tab_size*=this->dim(i);



    return MappedDenseLattice<data_type>(this->raw_data_ptr(),row_size, column_size, tab_size);






}

template<class T, size_t _order>
template<typename otherDerived, typename Op,typename index_param_type,typename boost::enable_if< internal::is_SparseMIA<otherDerived>, int >::type>
void  DenseMIA<T,_order>::merge(const MIA<otherDerived>  & restrict_libmia b,const Op& op,const std::array<index_param_type,_order>& index_order)
{


    #ifdef LIBMIA_CHECK_DIMS
    this->check_merge_dims(b,index_order);
    #endif

    static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Must use an array convertable to index_type");
    auto & b_derived=b.final_derived();
    for(auto it=b_derived.storage_begin();it<b_derived.storage_end();++it){
        //the index values of b may correspond to a shuffled version of the default linear index
        auto default_order_idx=b_derived.convert_to_default_linIdxSequence(b_derived.index_val(*it));
        //calculate the lhs_index based on how the MIA indices matched up
        auto lhs_index=internal::sub2ind_reorder(b_derived.ind2sub(default_order_idx),index_order,b_derived.dims());
        this->atIdx(lhs_index)=op(this->atIdx(lhs_index),this->convert(b_derived.data_val(*it)));

    }


}

template<class T, size_t _order>
template<typename otherDerived, typename Op,typename boost::enable_if< internal::is_SparseMIA<otherDerived>, int >::type>
void  DenseMIA<T,_order>::merge(const MIA<otherDerived>  & restrict_libmia b,const Op& op)
{


    #ifdef LIBMIA_CHECK_DIMS
    this->check_merge_dims(b);
    #endif

    auto & b_derived=b.final_derived();
    for(auto it=b_derived.storage_begin();it<b_derived.storage_end();++it){
        //the index values of b may correspond to a shuffled version of the default linear index
        auto default_order_idx=b_derived.convert_to_default_linIdxSequence(b_derived.index_val(*it));
        //calculate the lhs_index based on how the MIA indices matched up

        this->atIdx(default_order_idx)=op(this->atIdx(default_order_idx),this->convert(b_derived.data_val(*it)));

    }


}


/*! @} */

}

#endif // DENSEMIA_H_INCLUDED
