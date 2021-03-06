// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef MIA_H
#define MIA_H

#include <array>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "libdivide.h"

#include "Index.h"

#include "MIA_Expr.h"
#include "LibMIAUtil.h"
#include "IndexUtil.h"
#include "FunctionUtil.h"
namespace LibMIA
{

namespace internal
{
template<class Derived>
struct data_type<MIA<Derived> >: public data_type<Derived> {};

template<class Derived>
struct data_type_ref<MIA<Derived> >: public data_type_ref<Derived> {};

template<class Derived>
struct const_data_type_ref<MIA<Derived> >: public const_data_type_ref<Derived> {};

template<class Derived>
struct index_type<MIA<Derived> >: public index_type<Derived> {};

template<class Derived>
struct order<MIA<Derived> >: public order<Derived> {};

template<class Derived>
struct data_iterator<MIA<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct storage_iterator<MIA<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<MIA<Derived> >: public const_storage_iterator<Derived> {};

template<class Derived>
struct FinalDerived<MIA<Derived> >:public FinalDerived<Derived>{};

}


template<class MIA_Type,class UnaryOperation,
    typename  boost::enable_if<
        internal::is_MIA<
            typename std::remove_const<MIA_Type>::type
        >,int
    >::type=0
>
typename MIANonlinearFuncType<MIA_Type>::type do_func(const MIA_Type& mia, UnaryOperation op){
    typename MIANonlinearFuncType<MIA_Type>::type ret=mia;
    ret.apply_func(op);
    return ret;
}

template<class MIA_Type,
    typename  boost::enable_if<
        internal::is_MIA<
            typename std::remove_const<MIA_Type>::type
        >,int
    >::type=0
>
typename MIANonlinearFuncType<MIA_Type>::type sqrt(const MIA_Type& mia){
    typedef typename internal::data_type<MIA_Type>::type data_type;
    data_type (*fpc)(data_type) = &std::sqrt;
    return do_func(mia,std::ptr_fun(fpc));
}

//!raise MIA to a power
template<class MIA_Type,
    typename  boost::enable_if<
        internal::is_MIA<
            typename std::remove_const<MIA_Type>::type
        >,int
    >::type=0
>
typename MIANonlinearFuncType<MIA_Type>::type pow(const MIA_Type& mia,double exp){
    typedef typename internal::data_type<MIA_Type>::type data_type;
    data_type (*fpc)(data_type,data_type) = &std::pow;
    auto f1 = std::bind(fpc, std::placeholders::_1, exp);

    return do_func(mia,f1);
}

//!negate every non-zerp
template<class MIA_Type,
    typename  boost::enable_if<
        internal::is_MIA<
            typename std::remove_const<MIA_Type>::type
        >,int
    >::type=0
>
typename MIANonlinearFuncType<MIA_Type>::type negate(const MIA_Type& mia){
    typedef typename internal::data_type<MIA_Type>::type data_type;

    return do_func(mia,std::negate<data_type>());
}

//!Helper structure to ensure that anytime dimensions are changed, the dimensionality is automatically recalculated
template<typename index_type,size_t order>
struct dimsContainer{
    std::array<index_type,order> m_dims;
	std::array<size_t, order> temp_m_dims;
    index_type m_dimensionality;
    dimsContainer(){
        m_dims.fill(0);
        m_dimensionality=0;
    }
    void calculate_dimensionality(){
        m_dimensionality=1;
        for(auto i=this->m_dims.begin();i<this->m_dims.end();i++){
            m_dimensionality*=*i;

        }
    }

    dimsContainer(const std::array<index_type,order> & _dims){
        m_dims=_dims;
        calculate_dimensionality();

    }


    template<typename ...E>
    dimsContainer(E...e):temp_m_dims{{e...}}{
		std::copy(temp_m_dims.begin(), temp_m_dims.end(), m_dims.begin());
        calculate_dimensionality();
    }
    bool operator==(const std::array<index_type,order> & _dims)const{
        return m_dims==_dims;
    }
    bool operator!=(const std::array<index_type,order> & _dims)const{
        return !(m_dims==_dims);
    }
    index_type& operator[](int idx) {
        return m_dims[idx];
    }
    const index_type& operator[](int idx) const {
        return m_dims[idx];
    }
    const std::array<index_type,order> &  dims() const{
        return m_dims;
    }
    index_type dimensionality() const{
        return m_dimensionality;
    }
};

/** \addtogroup mia Multi-Index Array Classes
*  @{
*/

//!  Base class for multi-index array classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  multi-index array types. Provides common operations and functions.

  \tparam Derived   should only be DenseMIABase or SparseMIABase type.
*/
template
<

    class Derived
>
class MIA
{

public:

    typedef typename internal::index_type<MIA>::type index_type;
    typedef typename internal::data_type<MIA>::type data_type;
    typedef typename internal::data_type_ref<MIA>::type data_type_ref;
    typedef typename internal::const_data_type_ref<MIA>::type const_data_type_ref;
    typedef typename internal::FinalDerived<MIA>::type FinalDerived;
    typedef typename std::make_unsigned<index_type>::type unsigned_index_type;
    constexpr static size_t mOrder=internal::order<Derived>::value;
    typedef libdivide::divider<unsigned_index_type> fast_divisor;
    typedef std::array<unsigned_index_type,mOrder> accumulator_type;
    typedef std::array<unsigned_index_type,mOrder> multiplier_type;
    typedef std::array<fast_divisor,mOrder> fast_accumulator_type;

    Derived& derived() { return *static_cast<Derived*>(this); }
    /** \returns a const reference to the derived object */
    const Derived& derived() const { return *static_cast<const Derived*>(this); }

    FinalDerived& final_derived() {

        return derived().final_derived();
    }
    /** \returns a const reference to the derived object */
    const FinalDerived& final_derived() const {

        return derived().final_derived();
    }

    template<class...Ts>
    auto operator()(Ts...ts)->decltype(perform_unary<FinalDerived,typename internal::Indicial_Sequence<Ts...>::sequence>(this->final_derived())) {
        static_assert(internal::check_mia_indexing<MIA,Ts...>::type::value,"Number of dimensions must be same as <order> and each given index must be created using the MIAIndex macro.");
        return perform_unary<FinalDerived,typename internal::Indicial_Sequence<Ts...>::sequence>(final_derived());

        //return MIA_Atom<FinalDerived,typename internal::Indicial_Sequence<Ts...>::sequence,false>(&final_derived());

    }

    template<class...Ts>
    auto operator()(Ts...ts)const->decltype(perform_unary<const FinalDerived,typename internal::Indicial_Sequence<Ts...>::sequence>(this->final_derived())) {
        static_assert(internal::check_mia_indexing<MIA,Ts...>::type::value,"Number of dimensions must be same as <order> and each given index must be created using the MIAIndex macro.");
        return perform_unary<const FinalDerived,typename internal::Indicial_Sequence<Ts...>::sequence>(final_derived());

        //return MIA_Atom<FinalDerived,typename internal::Indicial_Sequence<Ts...>::sequence,false>(&final_derived());

    }

    MIA():m_dims(){

    }


    template<typename index_param_type>
    MIA(const std::array<index_param_type,mOrder > &_dims, SolveInfo _solveInfo=NoInfo):m_dims(_dims),mSolveInfo(_solveInfo){
        static_assert(internal::check_index_compatibility<index_type,index_param_type>::type::value,"Dimensions must be given in a data type convertable to index_type");


    }

    MIA(const std::array<index_type,mOrder > &_dims, SolveInfo _solveInfo=NoInfo): m_dims(_dims),mSolveInfo(_solveInfo) {}

    template<typename... Dims>
    MIA(Dims... dims):m_dims(dims...) {
        static_assert(internal::check_mia_constructor<MIA,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
    }








    void init(const std::array<index_type,mOrder> _dims,SolveInfo _solveInfo=NoInfo){
        m_dims(_dims);
        mSolveInfo(_solveInfo);
    }

//    template<class otherDerived,class index_param_type>
//    void assign(const MIA<otherDerived>& otherMIA,const std::array<index_param_type,mOrder>& index_order){
//        std::cout << "MIA Assign " << std::endl;
//        final_derived().assign(otherMIA.final_derived(),index_order);
//    }

    template<class otherDerived>
    FinalDerived& operator=(const MIA<otherDerived>& otherMIA){

        return final_derived()=otherMIA.final_derived();
    }

    //!Assignment operator. Will call Derived's operator
    FinalDerived& operator=(const MIA& otherMIA){

        return final_derived()=otherMIA.final_derived();
    }

    //!Assignment move operator. Will call Derived's operator
    FinalDerived& operator=(MIA&& otherMIA){

        return final_derived()=std::move(otherMIA.final_derived());
    }

    //!Applies the given function, component-wise to each data element. Typically some sort of nonlinear function. Can be lambda, std::function, or custom functor
    template<class UnaryOperation>
    void apply_func(UnaryOperation op){
        for(auto it=derived().data_begin(); it< derived().data_end();++it)
            *it=this->convert(op(*it));

    }

    void negate(){
        apply_func(std::negate<data_type>());
    }


//    toLatticeExpression(std::array<size_t> outer_product_indices, std::array<size_t> inner_product_indices,,std::array<size_t> inter_product_indices){
//        return toLatticeCopy(outer_product_indices, inner_product_indices, inter_product_indices );
//
//    }
//
    /** Sets all mia data to one.*/
    void ones(){
        std::fill ( derived().data_begin(), derived().data_end(), 1);
    }

    /** Sets all mia data to one.*/
    void zeros(){
        std::fill ( derived().data_begin(), derived().data_end(), 0);
    }

    /** Sets all mia data to given value.*/
    void fill(data_type val){
        std::fill ( derived().data_begin(), derived().data_end(), val);
    }

    index_type dim(size_t i) const{
        assert(i<mOrder);
        return m_dims[i];
    }

    //!  Sets MIA data to uniformly distributed random values.
    /*!
    If MIA is sparse, only elements already designated as nonzero are set to random values.
    Range is specified using parameters. Will throw a MIAParameterException exception if \f$low>high\f$.

    */
    void randu(int low=0, int high=1){
        using namespace boost::numeric;
        if (low>=high){
            throw MIAParameterException("Lower bound of random numbers must be stricly smaller than upper bound.");
        }
        boost::uniform_real<> uni_dist(low,high);
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(LibMIA_gen(), uni_dist);
        typedef boost::numeric::converter<data_type,boost::uniform_real<>::result_type> to_mdata_type;
        for (auto i=derived().data_begin();i<derived().data_end();++i){
            *i=to_mdata_type::convert(uni());
        }

    }

    //! Converts a scalar value to data_type
    /*!
        \tparam from_data_type the data_type you are converting from
    */
    template<class from_data_type>
    inline static data_type convert(const from_data_type from){

        return internal::convert<data_type,from_data_type>(from);
    }

    inline static data_type& convert(data_type& from){

        return from;
    }

    inline static const data_type& convert(const data_type& from){

        return from;
    }

    const std::array<index_type,internal::order<MIA>::value>& dims() const{
        return m_dims.dims();
    }
    index_type dimensionality() const{
        return m_dims.dimensionality();
    }
    void set_dims(const std::array<index_type,mOrder>& _dims){
        m_dims=_dims;

    }

    index_type sub2ind(const std::array<index_type,mOrder> & indices) const{
        return internal::sub2ind(indices,this->dims());
    }

    std::array<index_type,mOrder> ind2sub(index_type idx) const{
        return internal::ind2sub(idx, this->dims());
    }

    //! Returns scalar data at given indices
    /*!
        \param[in] indices variadic parameter. Will assert a compile error if size of indices!=mOrder or if Indices datatype are not convertible to index_type
    */
    template<typename... Indices>
    const data_type_ref at(Indices... indices) const {
        static_assert(internal::check_mia_constructor<MIA,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        std::array<index_type,internal::order<MIA>::value> temp = {{indices...}};
        return at(temp);
    }

    //! Returns scalar data at given indices
    /*!
        \param[in] indices variadic parameter. Will assert a compile error if size of indices!=mOrder or if Indices datatype are not convertible to index_type
    */
    template<typename... Indices>
    data_type_ref at(Indices... indices) {
        static_assert(internal::check_mia_constructor<MIA,Indices...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
        //auto temp = {{indices...}};
        std::array<index_type,internal::order<MIA>::value> temp = {{indices...}};
        return at(temp);
    }


    //! Returns scalar data at given indices
    const_data_type_ref at(const std::array<index_type, mOrder>& indices) const{

        return derived().atIdx(this->sub2ind(indices));

    }

    //! Returns scalar data at given indices
    data_type_ref at(const std::array<index_type, mOrder> & indices) {

        return derived().atIdx(this->sub2ind(indices));
    }

    //! Returns scalar data at given indices
    inline const_data_type_ref atIdx(index_type idx) const{
#ifdef LIBMIA_CHECK_DIMS
        if(idx<0||idx>=this->dimensionality())
            throw MIAParameterException("Index is out of range of the MIA.");
#endif
        return derived().atIdx(idx);

    }

    //! Returns scalar data at given indices
    inline data_type_ref atIdx(index_type idx) {
#ifdef LIBMIA_CHECK_DIMS
        if(idx<0||idx>=this->dimensionality())
            throw MIAParameterException("Index is out of range of the MIA.");
#endif
        return derived().atIdx(idx);
    }


    //!returns the data value at the internal DOK data container index. Same as atIdx for DenseMIAs, but for SparseMIAs, dok_index indexes only the nonzeros
    data_type_ref data_at(size_t dok_index)
    {
        return *(final_derived().data_begin()+dok_index);

    }
    //!returns the data value at the internal DOK data container index. Same as atIdx for DenseMIAs, but for SparseMIAs, dok_index indexes only the nonzeros
    const_data_type_ref data_at(size_t dok_index) const
    {
        return *(final_derived().data_begin()+dok_index);

    }




    template<size_t no_indices, size_t no_partitions>
    typename MIAUnaryType<Derived,no_indices>::type contract(const std::array<int,no_indices> & contract_indices,const std::array<int,no_partitions> & contract_partitions){

        return derived().contract_attract(contract_indices,contract_partitions,std::array<int,0>(),std::array<int,0>());
    }

//    template<size_t no_indices>
//    typename MIAUnaryType<Derived,no_indices>::type attract(const std::array<int,no_indices> & attract_indices){
//
//        return derived().contract_attract(std::array<int,0>(),attract_indices);
//    }

    SolveInfo solveInfo() const{
        return mSolveInfo;
    }
    void setSolveInfo(SolveInfo _solveInfo){
        mSolveInfo=_solveInfo;
    }
private:
    dimsContainer<index_type,mOrder> m_dims;


protected:










    SolveInfo mSolveInfo=NoInfo;



    template<class otherDerived,class index_param_type>
    void check_merge_dims(const MIA<otherDerived> &b,const std::array<index_param_type,mOrder>& index_order) const
    {

#ifdef LIBMIA_CHECK_DIMS
		for(size_t i=0;i<index_order.size();++i)
            if(b.dim(index_order[i])!=m_dims[i])
                throw MIParameterException("MIA dimensions must be identical for merger operation (+,-, etc).");
#endif

    }

    template<class otherDerived>
    void check_merge_dims(const MIA<otherDerived> &b) const
    {
#ifdef LIBMIA_CHECK_DIMS
        if(this->dims()!=b.dims())
            throw MIAParameterException("MIA dimensions must be identical for merger operation (+,-, etc).");
#endif

    }


};



/*! @} */


} //namespace LibMIA


#endif // MIA_H
