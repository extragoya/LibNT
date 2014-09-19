// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


// Implementation of the DenseLattice class. Underlying data
// structure and algorithms are those provided by the
// Eigen C++ Library.

#ifndef DENSELATTICEBASE_H
#define DENSELATTICEBASE_H
#include <iostream>
#include <fstream>
#include <string>


#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/type_traits/is_floating_point.hpp>
//#include <boost/timer/timer.hpp>

#include <Eigen/Dense>

#include "LibMIAUtil.h"
#include "Lattice.h"
#include "Packer.h"



namespace LibMIA
{

namespace internal
{



template<class Derived>
struct data_type<DenseLatticeBase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<DenseLatticeBase<Derived> >: public index_type<Derived> {};


template<class Derived>
struct data_iterator<DenseLatticeBase<Derived> >: public data_iterator<Derived> {};

template<class Derived>
struct const_data_iterator<DenseLatticeBase<Derived> >: public const_data_iterator<Derived> {};

template<class Derived>
struct storage_iterator<DenseLatticeBase<Derived> >: public storage_iterator<Derived> {};

template<class Derived>
struct const_storage_iterator<DenseLatticeBase<Derived> >: public const_storage_iterator<Derived> {};

template<class Derived>
struct Data<DenseLatticeBase<Derived> >: public Data<Derived> {};

}

/** \addtogroup lattice Lattice Classes
 *  @{
 */

//!  Base class for dense lattice classes.
/*!
  This class acts as the base class for the parametric subclass pattern,
  which is more commonly known as the CTRP. It is the base class for all
  dense lattice types. Provides operations and functions common to all dense
  lattices.

  \tparam Derived   should only be a dense lattice class.
*/
template <class Derived>
class DenseLatticeBase: public Lattice<DenseLatticeBase<Derived > >//, boost::multipliable<DenseLattice<T> >
{
public:


    typedef typename internal::data_type<DenseLatticeBase>::type data_type;
    typedef typename internal::index_type<DenseLatticeBase>::type index_type;
    typedef typename internal::Data<DenseLatticeBase>::type Data;
    typedef typename internal::data_iterator<DenseLatticeBase>::type data_iterator;
    typedef typename internal::const_data_iterator<DenseLatticeBase>::type const_data_iterator;
    typedef typename internal::storage_iterator<DenseLatticeBase>::type storage_iterator;
    typedef typename internal::const_storage_iterator<DenseLatticeBase>::type const_storage_iterator;
    typedef typename Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> > matrix_type;
    typedef typename Eigen::Map<const Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> > const_matrix_type;
    typedef typename Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, 1> > vector_type;
    typedef typename Eigen::Map<const Eigen::Matrix<data_type, Eigen::Dynamic, 1> > const_vector_type;

    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }
    /** \returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

    DenseLatticeBase(){}

    //!  Performs lattice product.
    /*!
        c=a*b. Performs matrix/matrix product across each of the operand types. Uses <a href="http://eigen.tuxfamily.org">Eigen</a> to perform
        matrix products. Operation is only defined for two dense lattice of the same scalar type (which may change as EigenLib evolves).

        \returns DenseLattice<data_type> regardless of the operand types.
    */
    template<class otherDerived>
    typename DenseProductReturnType<Derived,otherDerived>::type
    operator*(const DenseLatticeBase<otherDerived> &b) const;


    //!  Performs lattice product with sparse lattice.
    /*!
        c=a*b. Performs matrix/matrix product across each of the operand types. Operation is only defined for two dense lattice of the same scalar type (which may change as EigenLib evolves).

        \returns DenseLattice<data_type> regardless of the operand types.
    */
    template<class otherDerived>
    typename SparseProductReturnType<Derived,otherDerived>::type
    operator*(SparseLatticeBase<otherDerived> &b) const;


    //!  Solves linear system of equations.
    /*!
    x=a.solve(b), solves a sequence of systems of linear equations represented by a*x=b. Uses <a href="http://eigen.tuxfamily.org">Eigen</a> to compute
        the solution. System can be over-determined. Operation is only defined for two dense lattice of the same scalar type (which may change as EigenLib evolves). Scalar type
        must be floating point.

        \returns DenseLattice<data_type> regardless of the operand types.
    */
    template<class otherDerived>
    typename DenseSolveReturnType<Derived,otherDerived>::type
    solve(const DenseLatticeBase<otherDerived> &b) const;


    //!  Solves linear system of equations.
    /*!
    x=a.solve(b), solves a sequence of systems of linear equations represented by a*x=b, where b is sparse. Uses <a href="http://eigen.tuxfamily.org">Eigen</a> to compute
        the solution. System can be over-determined. Operation is only defined for two lattices of the same scalar type (which may change as EigenLib evolves). Scalar type
        must be floating point.

        \returns DenseLattice<data_type> regardless of the operand types.
    */
    template<class otherDerived>
    typename DenseSolveReturnType<Derived,otherDerived>::type
    solve(SparseLatticeBase<otherDerived> &b) const;



    void print() const;

    bool save(const std::string & _filename) const;

    DenseLattice<data_type> transpose() const;
    void inPlaceTranspose();

    template<class otherDerived>
    bool operator==(const DenseLatticeBase<otherDerived> &b) const;

    template<class otherDerived>
    bool operator==(SparseLatticeBase<otherDerived> &b){
        return b==*this;
    }

    template<class otherDerived>
    bool fuzzy_equals(SparseLatticeBase<otherDerived> &b,double precision){
        return b.fuzzy_equals(*this,precision);
    }

    template<class otherDerived>
    bool fuzzy_equals(const DenseLatticeBase<otherDerived> &b,double precision) const;

    template<class otherDerived>
    bool operator!=(SparseLatticeBase<otherDerived> &b){
        return !(b==*this);
    }



    //!  Sets each tab to an identity matrix.
    /*!
        If non-square, it follows <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1MatrixBase.html#a0650b65c6ae6c3d19a138b72a6d68568">this</a> format.
    */
    void eye();

    //!  Destructive add.
    /*!
        \returns DenseLattice<data_type>. If data_type is a subtype of b's data_type, precision may be lost. Conversion executed using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    DenseLatticeBase & operator+=(const DenseLatticeBase<otherDerived> &b)
    {
        this->check_merge_dims(b);
        std::plus<data_type> op;
        merge(b,op);
        return *this;

    }

    //!  Destructive subtract.
    /*!
        \returns DenseLattice<data_type>. If data_type is a subtype of b's data_type, precision may be lost. Conversion executed using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    DenseLatticeBase & operator-=(const DenseLatticeBase<otherDerived> &b)
    {
        this->check_merge_dims(b);
        std::minus<data_type> op;
        merge(b,op);
        return *this;

    }

    //!  Binary addition.
    /*!
        \returns DenseLattice<super_data_type>. Super_data_type is determined using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    typename DenseMergeReturnType<Derived,otherDerived>::type operator+(DenseLatticeBase<otherDerived> &b)
    {
        typedef typename DenseMergeReturnType<Derived,otherDerived>::type RType;
        this->check_merge_dims(b);
        RType c(*this);
        c+=b;
        return c;

    }

    //!  Binary subtraction.
    /*!
        \returns DenseLattice<super_data_type>. Super_data_type is determined using
        <a href=" http://www.boost.org/doc/libs/1_34_1/libs/numeric/conversion/doc/index.html">Boost's Conversion Library</a>.
    */
    template<class otherDerived>
    typename DenseMergeReturnType<Derived,otherDerived>::type operator-(DenseLatticeBase<otherDerived> &b)
    {
        typedef typename DenseMergeReturnType<Derived,otherDerived>::type RType;
        this->check_merge_dims(b);
        RType c(*this);
        c-=b;
        return c;

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

    storage_iterator begin()
    {


        return derived().begin();
    }

    storage_iterator end()
    {


        return derived().end();
    }
    const_storage_iterator begin() const
    {


        return derived().begin();
    }

    const_storage_iterator end() const
    {


        return derived().end();
    }

    //!  Wraps specified tab with an <a href="http://eigen.tuxfamily.org">Eigen</a> matrix.
    matrix_type tab_matrix(int _tab){
        return matrix_type(&(derived()(0,0,_tab)),this->height(),this->width());

    }

    const_matrix_type tab_matrix(int _tab) const{
        return const_matrix_type( const_cast<data_type *>(&(derived()(0,0,_tab))),this->height(),this->width());

    }

    //!  Wraps specified column with an <a href="http://eigen.tuxfamily.org">Eigen</a> vector.
    vector_type column_vector(int _column,int _tab){
        return vector_type(&(derived()(0,_column,_tab)),this->height());

    }

    const_vector_type column_vector(int _column,int _tab) const{
        return const_vector_type(&(derived()(0,_column,_tab)),this->height());

    }

    const data_type& operator()(index_type i, index_type j, index_type k) const
    {

        return atIdx(this->sub2ind(i,j,k));
    }
    data_type& operator()(index_type i, index_type j, index_type k)
    {

        return atIdx(this->sub2ind(i,j,k));
    }


    //! Returns scalar data at given linear index
    const data_type& atIdx(index_type idx) const{

        //return lin index
        return *(derived().data_begin()+idx);
    }

    //! Returns scalar data at given linear index
    data_type& atIdx(index_type idx) {

        //return lin index
        return *(derived().data_begin()+idx);
    }

protected:

    template<class otherDerived,class Predicate>
    bool compare_with_dense(const DenseLatticeBase<otherDerived>& b, Predicate pred) const;

    template<class otherDerived,bool LSQR>
    inline void perform_solve(SparseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const;

    template<class otherDerived,bool LSQR>
    inline void perform_solve(const DenseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const;



    template <class otherDerived,class Op>
    void merge(const DenseLatticeBase<otherDerived> &b,Op op)
    {


        typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;
        typename otherDerived::const_storage_iterator other_it=b.begin();
        for(storage_iterator i=begin(); i<end(); i++)
            *i=to_mdata_type::convert(op(*i,*other_it++));

    }
private:



};



template <class Derived>
template<class otherDerived>
inline
typename DenseProductReturnType<Derived,otherDerived>::type
DenseLatticeBase<Derived>::operator*(const DenseLatticeBase<otherDerived> & restrict_libmia b) const restrict_libmia
{


    this->check_mult_dims(b);


    typedef typename DenseProductReturnType<Derived,otherDerived>::type RType;


    RType c(this->height(),b.width(),this->depth());

    for (int i=0; i<this->depth(); i++)
    {

        c.derived().tab_matrix(i)=(this->derived().tab_matrix(i))*(b.derived().tab_matrix(i));


    }

    return c;
}


template <class Derived>
template<class otherDerived>
inline
typename SparseProductReturnType<Derived,otherDerived>::type
DenseLatticeBase<Derived>::operator*(SparseLatticeBase<otherDerived> & restrict_libmia b) const restrict_libmia
{



    //boost::timer::cpu_timer timer,timer_minor;
    this->check_mult_dims(b);
    auto & a=*this;
    typedef typename SparseProductReturnType<Derived,otherDerived>::type c_type;
    typedef typename Derived::index_type a_index_type;

    typename c_type::Indices c_indices;
    c_indices.reserve(b.size());
    typename c_type::Data c_data;
    c_data.reserve(b.size());
    typename internal::data_type<c_type>::type cur_c_data;
    //timer_minor=boost::timer::cpu_timer();

    b.sort();

    auto b_temp_begin=b.index_begin();
    auto b_temp_end=b_temp_begin;

    while(b_temp_begin<b.index_end()){
        auto cur_tab=b.tab(*b_temp_begin);
        auto cur_col=b.column(*b_temp_begin);
        b_temp_end=b_temp_begin+1;
        while(b_temp_end<b.index_end()&&b.tab(*b_temp_end)==cur_tab && b.column(*b_temp_end)==cur_col){
            b_temp_end++;
        }
        for(a_index_type a_rows=0;a_rows<a.height();++a_rows){
            cur_c_data=0;
            for(auto b_it=b_temp_begin;b_it<b_temp_end;++b_it){
                cur_c_data+=b.data_at(b_it-b.index_begin())*a(a_rows,b.row(*b_it),cur_tab);

            }
            if(std::abs(cur_c_data)>SparseTolerance<decltype(cur_c_data)>::tolerance){
                c_data.push_back(cur_c_data);
                c_indices.push_back(a_rows+(cur_col+cur_tab*b.width())*a.height());
            }
        }
        b_temp_begin=b_temp_end;
    }
    c_type ret(std::move(c_data),std::move(c_indices),a.height(),b.width(),this->depth());
    ret.set_linIdxSequence(ColumnMajor);
    //std::cout << "dense x sparse time " << boost::timer::format(timer.elapsed()) << std::endl;
    return ret;





}

template <class Derived>
auto DenseLatticeBase<Derived>::transpose()  const->DenseLattice<data_type>{



    DenseLattice<data_type> c(this->width(),this->height(),this->depth());


    for (int i=0; i<this->depth(); i++)
    {

        c.derived().tab_matrix(i)=tab_matrix(i).transpose();

    }

    return c;


}

template <class Derived>
void DenseLatticeBase<Derived>::inPlaceTranspose() {





    for (int i=0; i<this->depth(); i++)
    {
        this->derived().tab_matrix(i).transposeInPlace();

    }
    std::swap(this->m_height,this->m_width);



}

namespace {
//!Helper class to pull least squares or householder inversion
template<bool LSQR,class Lattice>
struct dense_lattice_solver;

template<class Lattice>
struct dense_lattice_solver<false,Lattice> {

    static auto get_solver(const Lattice & lattice,int _tab,SolveInfo &_solveInfo)->decltype(lattice.tab_matrix(_tab).fullPivLu()){
        auto _solver=lattice.tab_matrix(_tab).fullPivLu();
        if(!_solver.isInvertible())
            _solveInfo=RankDeficient;
        else
            _solveInfo=FullyRanked;
        return _solver;

    }

};

template<class Lattice>
struct dense_lattice_solver<true,Lattice> {

    static auto get_solver(const Lattice & lattice,int _tab,SolveInfo &_solveInfo)->decltype(lattice.tab_matrix(_tab).jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)){
        auto _SVD=lattice.tab_matrix(_tab).jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        if(_SVD.nonzeroSingularValues()!=lattice.width())
            _solveInfo=RankDeficient;
        else
            _solveInfo=LeastSquares;
        return _SVD;
    }

};

}

template <class Derived>
template<class otherDerived>
inline
typename DenseSolveReturnType<Derived,otherDerived>::type
DenseLatticeBase<Derived>::solve(const DenseLatticeBase<otherDerived> &b) const{


    this->check_solve_dims(b);



    typedef typename DenseSolveReturnType<Derived,otherDerived>::type RType;

    RType c(this->width(),b.width(),this->depth());

    if (this->width()==this->height()){
        perform_solve<otherDerived,false>(b,c);
    }
    else if (this->width()<this->height()){
       perform_solve<otherDerived,true>(b,c);

    }



    return c;


}

template <class Derived>
template<class otherDerived>
inline
typename DenseSolveReturnType<Derived,otherDerived>::type
DenseLatticeBase<Derived>::solve(SparseLatticeBase<otherDerived> &b) const{


    this->check_solve_dims(b);



    typedef typename DenseSolveReturnType<Derived,otherDerived>::type RType;

    RType c(this->width(),b.width(),this->depth());

    if (this->width()==this->height()){
        perform_solve<otherDerived,false>(b,c);
    }
    else if (this->width()<this->height()){
        perform_solve<otherDerived,true>(b,c);

    }



    return c;


}


template <class Derived>
template<class otherDerived,bool LSQR>
inline void DenseLatticeBase<Derived>::perform_solve(SparseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const{



    typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> temp_vector_type;

    b.sort();


    SolveInfo _solveInfo;
    c.setSolveInfo(FullyRanked);
    auto b_temp_tab_begin=b.index_begin();
    decltype(b_temp_tab_begin) b_temp_tab_end;
    auto b_index_end=b.index_end();
    temp_vector_type temp_column(temp_vector_type::Zero(b.height()));
    for (index_type cur_tab=0; cur_tab<this->depth(); ++cur_tab)
    {


        b_temp_tab_begin=b.find_tab_start_idx(cur_tab,b_temp_tab_begin,b_index_end);
        if(b_temp_tab_begin==b_index_end)
            break;

        b_temp_tab_end=b.find_tab_end_idx(cur_tab,b_temp_tab_begin,b_index_end);


        if(b_temp_tab_end==b_temp_tab_begin)
            continue;

        auto _solver=dense_lattice_solver<LSQR,DenseLatticeBase>::get_solver(*this,cur_tab,_solveInfo);
        c.setSolveInfo(_solveInfo);

        auto b_temp_col_begin=b_temp_tab_begin;
        while(b_temp_col_begin<b_temp_tab_end){
            auto cur_col=b.column(*b_temp_col_begin);
            auto b_temp_col_end=b_temp_col_begin;
            temp_column.setZero();
            while(b_temp_col_end < b_temp_tab_end && b.column(*b_temp_col_end)==cur_col){
                temp_column(b.row(*b_temp_col_end))=b.data_at(b_temp_col_end-b.index_begin()); //convert to rows
                ++b_temp_col_end;
            }
            c.derived().tab_matrix(cur_tab).col(cur_col)=_solver.solve(temp_column);

            b_temp_col_begin=b_temp_col_end;
        }

        b_temp_tab_begin=b_temp_tab_end;

    }



}




template <class Derived>
template<class otherDerived,bool LSQR>
inline void DenseLatticeBase<Derived>::perform_solve(const DenseLatticeBase<otherDerived> &b, typename DenseSolveReturnType<Derived,otherDerived>::type & c) const{




    c.setSolveInfo(NoInfo);
    SolveInfo _solveInfo=NoInfo;
    for (int i=0; i<this->depth(); i++)
    {

        auto _solver=dense_lattice_solver<LSQR,DenseLatticeBase>::get_solver(*this,i,_solveInfo);
		if (_solveInfo==LibMIA::RankDeficient) //ensures that if we get only one tab Rank Deficient, we'll still set the entire solution as Rank Deficient
			c.setSolveInfo(LibMIA::RankDeficient);
        c.derived().tab_matrix(i)=_solver.solve(b.derived().tab_matrix(i));

    }
	if (c.solveInfo()!=LibMIA::RankDeficient) //if no tabs were RankDeficient, then set solveInfo to be either FullyRanked or LeastSquares
		c.setSolveInfo(_solveInfo);

}




template<class Derived>
void DenseLatticeBase<Derived>::print() const
{
    for (int i=0; i<this->depth(); i++)
    {


        std::cout<< "Tab: " << i << "\n";
        std::cout<< tab_matrix(i)<< "\n";


    }

}

template<class Derived>
void DenseLatticeBase<Derived>::eye()
{
    for (int i=0; i<this->depth(); i++)
    {


        tab_matrix(i).setIdentity();


    }

}

template<class Derived>
template<class otherDerived,class Predicate>
bool DenseLatticeBase<Derived>::compare_with_dense(const DenseLatticeBase<otherDerived> &b,Predicate predicate) const{

    if(this->height()!=b.height()||this->width()!=b.width()||this->depth()!=b.depth())
        return false;

    for(auto it=derived().data_begin(), itb=b.data_begin();it<derived().data_end();++it,++itb)
        if(!predicate(*it,*itb))
            return false;

    return true;


}

template<class Derived>
template<class otherDerived>
bool DenseLatticeBase<Derived>::operator==(const DenseLatticeBase<otherDerived> &b) const{

    typedef typename internal::index_type<otherDerived>::type b_index_type;
    std::function<bool(index_type,b_index_type)>pred=[](index_type lhs, b_index_type rhs){
        return lhs==rhs;
    };
    return compare_with_dense(b,pred);

}

template<typename Derived>
template<typename otherDerived>
bool DenseLatticeBase<Derived>::fuzzy_equals(const DenseLatticeBase<otherDerived> & b,double precision ) const
{


    typedef typename otherDerived::index_type b_index_type;
    std::function<bool(index_type,b_index_type)>pred=[precision](index_type lhs, b_index_type rhs){
        return isEqualFuzzy(lhs,rhs,precision);
    };
    return compare_with_dense(b,pred);

}


template<class Derived>
bool DenseLatticeBase<Derived>::save(const std::string & _filename) const
{
    std::ofstream outfile (_filename.c_str(),std::ios_base::binary);
    if (outfile.fail())
        return false;
    uint64_t temp;
    temp=pack754_64((pack754_type)this->height());
    outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));
    temp=pack754_64((pack754_type)this->width());
    outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));
    temp=pack754_64((pack754_type)this->depth());
    outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));

    for(auto it=derived().data_begin();it<derived().data_end();++it){
        temp=pack754_64(*it);
        outfile.write(reinterpret_cast<char *>(&temp),sizeof(uint64_t));
    }

    outfile.close();

    return true;

}




/*! @} */

}



#endif // DENSELATTICEBASE_H
