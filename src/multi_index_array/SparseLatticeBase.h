
#ifndef SPARSELATTICEBASE_H
#define SPARSELATTICEBASE_H
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/numeric/conversion/converter.hpp>



#include <Eigen/Sparse>
#include <Eigen/SuperLUSupport>
//

#include "Lattice.h"
#include "Util.h"
#include "tupleit.hh"

namespace LibMIA
{

namespace internal
{


template<class Derived>
struct data_type<SparseLatticeBase<Derived> >: public data_type<Derived> {};

template<class Derived>
struct index_type<SparseLatticeBase<Derived> >: public index_type<Derived> {};


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
        for (storage_iterator i=begin(); i<end(); i++)
        {

            full_tuple temp=*i;
            //T temp2=boost::get<0>(temp);

            std::cout << boost::get<0>(temp) <<"\t" ;
            std::cout << row(index(temp)) << "\t"<< column(index(temp)) << "\t" << tab(index(temp)) <<"\n";
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




        const index_type& (SparseLatticeBase::*index)(const_full_tuple) const = &SparseLatticeBase::index;
        storage_iterator temp=begin();

        if (m_sort_order!=sort_order || !is_sorted())
        {

            if (sort_order==ColumnMajor)
                std::sort(begin(),end(),boost::bind(&SparseLatticeBase::idx_less, this,boost::bind(index,this,_1),boost::bind(index,this,_2)));

            else
                std::sort(begin(),end(),boost::bind(&SparseLatticeBase::idx_less_rowmajor, this,boost::bind(index,this,_1),boost::bind(index,this,_2)));

            m_is_sorted=true;
            m_sort_order=sort_order;
        }


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


    index_type& index(full_tuple a)
    {
        return boost::get<1>(a);

    }

    const index_type& index(const_full_tuple a) const
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

    const data_type& data_at(index_type nnz_index) const
    {
        return *(data_begin()+nnz_index);

    }

    data_type& data_at(index_type nnz_index)
    {
        return *(data_begin()+nnz_index);

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



    storage_iterator begin();

    const_storage_iterator begin() const;

    storage_iterator end();

    const_storage_iterator end() const;


    //!  Utility function for mapping tabs to column-major matrices to use in lattice mulitplication
    /*!
    \param[in] row_map A sorted vector mapping compressed indices to full row indices of the current tab
    \param[in] inner_indices A sorted vector mapping shared column indices to full column indices of the current tab
    \param[in] t_begin t_end Storage iterator to the nonzeros of the current tab
    \param[out] mat Resulting CSC matrix. The CSC matrix should already have allocated vectors for data and indices.
    */
    template<class super_data_type>
    void to_matrix(std::vector<index_type> &row_map,std::vector<index_type> &inner_indices,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_cm& mat) const;


    //!  Utility function for mapping tabs to row-major matrices to use in lattice mulitplication. See to_matrix for more information.
    template<class super_data_type>
    void to_matrix_rowmajor(std::vector<index_type> &inner_indices,std::vector<index_type> &column_map,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_rm& mat) const;

private:


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
SparseLatticeBase<Derived>::begin()
{


    return iterators::makeTupleIterator(data_begin(),index_begin());



}


template<typename Derived>
inline typename SparseLatticeBase<Derived>::const_storage_iterator
SparseLatticeBase<Derived>::begin() const
{


    const_data_iterator start1=data_begin();
    const_index_iterator start2=index_begin();
    return boost::make_zip_iterator(boost::make_tuple(start1, start2));



}


template<typename Derived>
inline typename SparseLatticeBase<Derived>::storage_iterator
SparseLatticeBase<Derived>::end()
{


    return iterators::makeTupleIterator(data_end(),index_end());



}


template<typename Derived>
inline typename SparseLatticeBase<Derived>::const_storage_iterator
SparseLatticeBase<Derived>::end() const
{


    const_data_iterator start1=data_end();
    const_index_iterator start2=index_end();
    return boost::make_zip_iterator(boost::make_tuple(start1, start2));



}



//
//
template <class Derived>
template <class otherDerived>
typename SparseProductReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::operator*(SparseLatticeBase<otherDerived> &b)
{

    check_mult_dims(b);


    typedef typename ScalarPromoteType<Derived,otherDerived>::type super_data_type;

    typedef typename SparseProductReturnType<Derived,otherDerived>::type RType;
    typedef typename internal::index_type<otherDerived>::type b_index_type; //should be the same as index type, but just in case that changes in future versions

    otherDerived b_d=b.derived();
    Derived a_d=this->derived();
    sort(ColumnMajor); //tab/column major for A
    bool old_sort_order=b.sort_order();
    b.sort(RowMajor); //tab/row major for B

    //iterators for indices and data
    storage_iterator a_begin=begin();
    index_iterator a_index_begin=a_d.index_begin();
    data_iterator a_data_begin=data_begin();

    auto b_begin=b.begin();
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
        a_temp_end=std::upper_bound(a_temp_begin,a_index_end,k,boost::bind(std::less<index_type>(), _1, boost::bind(tab,this,_2)));
        b_temp_end=std::upper_bound(b_temp_begin,b_index_end,k,boost::bind(std::less<b_index_type>(), _1, boost::bind(tab_other,&b_d,_2)));

        //if both tabs have nonzeros, then perform matrix multiplication
        if (a_temp_end!=a_temp_begin && b_temp_end!=b_temp_begin)
        {

            //get unique columns of A //***todo - just use a binary predicate in std::unique_copy, avoiding having to store needless values
            a_columns.resize(a_temp_end-a_temp_begin);
            std::transform(a_temp_begin,a_temp_end,a_columns.begin(),boost::bind(column, this,_1));
            a_column_end=std::unique(a_columns.begin(),a_columns.end());
            a_columns.resize( a_column_end - a_columns.begin() );


            //get unique and sorted rows of A. This creates a map from compressed indices to actual ones //***todo - just a series of set unions to create a unique row index
            a_rows.resize(a_temp_end-a_temp_begin);
            std::transform(a_temp_begin,a_temp_end,a_rows.begin(),boost::bind(row, this,_1));
            std::sort(a_rows.begin(),a_rows.end());
            a_row_end=std::unique(a_rows.begin(),a_rows.end());
            a_rows.resize( a_row_end - a_rows.begin() );



            //get unique rows of B and a map from compressed indices to actual ones
            b_rows.resize(b_temp_end-b_temp_begin);
            std::transform(b_temp_begin,b_temp_end,b_rows.begin(),boost::bind(other_row, &b_d,_1));
            b_row_end=std::unique(b_rows.begin(),b_rows.end());
            b_rows.resize(b_row_end - b_rows.begin());

            //get unique and sorted columns of B. This creates a map from compressed indices to actual ones
            b_columns.resize(b_temp_end-b_temp_begin);
            std::transform(b_temp_begin,b_temp_end,b_columns.begin(),boost::bind(other_column, &b_d,_1));
            std::sort(b_columns.begin(),b_columns.end());
            b_column_end=std::unique(b_columns.begin(),b_columns.end());
            b_columns.resize( b_column_end - b_columns.begin() );


            //find union of A columns and B rows
            inner_indices.resize(a_columns.size()+b_rows.size());
            inner_end=std::set_union(a_columns.begin(),a_columns.end(),b_rows.begin(),b_rows.end(),inner_indices.begin());
            inner_indices.resize(inner_end-inner_indices.begin());



            //create temporary sparse matrices for current tab and calculate matrix product
            a_columns.resize(inner_indices.size());
            std::vector<index_type> rows_place_holder;
            rows_place_holder.resize(a_temp_end-a_temp_begin);
            MappedSparseMatrix_cm A(a_rows.size(),inner_indices.size(),a_temp_end-a_temp_begin,&a_columns[0],&rows_place_holder[0],&(*(a_data_begin()+(a_temp_begin-a_index_begin))));
            to_matrix(a_rows,inner_indices,a_begin+(a_temp_begin-a_index_begin),a_begin+(a_temp_end-a_index_begin),A);

            b_rows.resize(inner_indices.size());
            std::vector<b_index_type> columns_place_holder;
            columns_place_holder.resize(b_temp_end-b_temp_begin);
            MappedSparseMatrix_rm B(inner_indices.size(),b_columns.size(),b_temp_end-b_temp_begin,&b_rows[0],&columns_place_holder[0],&(*(b_data_begin()+(b_temp_begin-b_index_begin))));
            b.to_matrix_rowmajor(inner_indices,b_columns,b_begin+(b_temp_begin-b_index_begin),b_begin+(b_temp_end-b_index_begin),B);
            SparseMatrix_cm C=A*B;
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


        }
        a_temp_begin=a_temp_end;
        b_temp_begin=b_temp_end;

    }

    //resort B in the proper precedence
    b_d.sort(old_sort_order);
    //return sparse lattice using the data and indices vectors
    return RType(c_data,c_indices,this->height(),b.width(),this->depth());
}

template <class Derived>
template <class super_data_type>
inline void SparseLatticeBase<Derived>::to_matrix(std::vector<index_type> &row_map, std::vector<index_type> &inner_indices,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_cm& mat) const
{

    typedef boost::numeric::converter<super_data_type,data_type> toSuper;
    bool (SparseLatticeBase::*column_compare) (index_type, index_type) const = &SparseLatticeBase::column_compare;
    const index_type& (SparseLatticeBase::*index_bind)(const_full_tuple) const = &SparseLatticeBase::index;


    index_type outer_val=0;
    auto mat_rows=mat._innerIndexPtr();
    auto mat_columns=mat._outerIndexPtr();
    typename std::vector<index_type>::const_iterator cur_row;
    for(auto cur_it=inner_indices.begin(); cur_it<inner_indices.end(); cur_it++)
    {

        *mat_columns++=outer_val;
        cur_row=row_map.begin();
        while (column(index(*t_begin++))<*cur_it);


        while (column(index(*t_begin++))==*cur_it)
        {
            cur_row=std::lower_bound(cur_row,row_map.end(),row(index(*t_begin)));
            *mat_rows++=cur_row-row_map.begin();
            outer_val++;

        }



    }
    *mat_columns=outer_val;


}


template <class Derived>
template <class super_data_type>
inline void SparseLatticeBase<Derived>::to_matrix_rowmajor(std::vector<index_type> &inner_indices,std::vector<index_type> &column_map,storage_iterator t_begin, storage_iterator t_end, MappedSparseMatrix_rm& mat) const
{

    typedef boost::numeric::converter<super_data_type,data_type> toSuper;
    bool (SparseLatticeBase::*row_compare) (index_type, index_type) const = &SparseLatticeBase::row_compare;
    const index_type& (SparseLatticeBase::*index_bind)(const_full_tuple) const = &SparseLatticeBase::index;

    index_type outer_val=0;
    auto mat_columns=mat._innerIndexPtr();
    auto mat_rows=mat._outerIndexPtr();
    typename std::vector<index_type>::const_iterator cur_column;
    for(auto cur_it=inner_indices.begin(); cur_it<inner_indices.end(); cur_it++)
    {

        *mat_rows++=outer_val;
        cur_column=column_map.begin();
        while (row(index(*t_begin++))<*cur_it);

        while (row(index(*t_begin++))==*cur_it)
        {
            cur_column=std::lower_bound(cur_column,column_map.end(),column(index(*t_begin)));
            *mat_columns++=cur_column-column_map.begin();
            outer_val++;
        }


    }
    *mat_rows=outer_val;


}

template <class Derived>
template <class otherDerived>
typename SparseSolveReturnType<Derived,otherDerived>::type SparseLatticeBase<Derived>::solve(SparseLatticeBase<otherDerived> &b)
{

    check_solve_dims(b);



    typedef typename SparseSolveReturnType<Derived,otherDerived>::type c_type;
    typedef typename c_type::vector_type c_vector_type;

    typedef typename internal::index_type<otherDerived>::type b_index_type;
    typedef typename internal::storage_iterator<otherDerived>::type other_storage_iterator;
    typedef typename internal::index_iterator<otherDerived>::type other_index_iterator;
    //scalar datatype converter
    typedef boost::numeric::converter<data_type,typename internal::data_type<otherDerived>::type> to_mdata_type;

    otherDerived b_d=b.derived();
    Derived a_d=this->derived();
    sort(ColumnMajor); //tab/column major for A
    b.sort(ColumnMajor); //tab/column major for B

    //iterators for for indices and data
    storage_iterator a_begin=begin();
    index_iterator a_index_begin=index_begin();
    index_iterator a_index_end=index_end();
    index_iterator a_temp_begin=a_index_begin;
    index_iterator a_temp_end;

    other_storage_iterator b_begin=b.begin();
    other_index_iterator b_index_begin=b.index_begin();
    other_index_iterator b_temp_begin=b_index_begin;
    other_index_iterator b_index_end=b.index_end();
    other_index_iterator b_temp_end;


    //hold compressed indices for Sparse matrix of each tab
    std::vector<index_type> a_columns;
    a_columns.resize(this->width()+1);
    std::vector<index_type> a_rows;



    typedef typename Eigen::MappedSparseMatrix<data_type,Eigen::ColMajor,index_type> MappedSparseMatrix_cm; //Sparse Matrix type for A
    typedef Eigen::SparseLU<MappedSparseMatrix_cm,Eigen::SuperLU> LU_decomp;    //LU decomposition type for A


    c_type c(this->width(),b.width(),this->depth());   //create dense lattice to return and allocate memory


    typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> mapped_vector; //set to *this's datatype - conversion happens if needed
    mapped_vector b_vector(b.height());

    //create temporaries for SparseLattice member functions to use with boost::bind
    index_type (SparseLatticeBase::*column)(index_type lin_index) const = &SparseLatticeBase::column;
    index_type (otherDerived::*column_other)(index_type lin_index) const = &otherDerived::column;
    index_type (SparseLatticeBase::*row)(index_type lin_index) const = &SparseLatticeBase::row;
    index_type (SparseLatticeBase::*tab)(index_type lin_index) const = &SparseLatticeBase::tab;
    index_type (otherDerived::*tab_other)(index_type lin_index) const = &otherDerived::tab;
    const index_type& (otherDerived::*index_bind_other)(const_full_tuple) const = &otherDerived::index;


    for (int k=0; k<this->depth(); k++) //loop through every tab
    {


        //find the last occurence of current tab in both lattices
        a_temp_end=std::upper_bound(a_temp_begin,a_index_end,k,boost::bind(std::less<index_type>(), _1, boost::bind(tab,this,_2)));
        b_temp_end=std::upper_bound(b_temp_begin,b_index_end,k,boost::bind(std::less<index_type>(), _1, boost::bind(tab_other,&b_d,_2)));

        //tab from a must have nonzeros
        if (a_temp_end!=a_temp_begin)
        {

            //get rows of current tab of A
            a_rows.resize(a_temp_end-a_temp_begin);
            std::transform(a_temp_begin,a_temp_end,a_rows.begin(),boost::bind(row, this,_1));

            //created CCS matrix of current tab
            index_iterator cur_spot=a_temp_begin;
            //loop through every column of A
            for (int j=0; j<this->width(); j++)
            {
                //find the upper bound of values in the current column
                cur_spot=std::lower_bound(cur_spot,a_temp_end,j,boost::bind(std::less<index_type>(),boost::bind(column,this,_1), _2));
                if (cur_spot==a_index_end)
                {
                    std::stringstream t;
                    t << "Rank deficient tab. Column " << j << " in tab " << k << " of LHS has zero entries.";
                    throw RankDeficientException(t.str());

                }
                //record this location to use as outer index array for CCS
                a_columns[j] = cur_spot-a_temp_begin;


            }
            //end value of outer index array
            a_columns.back()=a_rows.size();

            //created a CCS matrix by mapping row and column vectors and also the pre-existing data of *this lattice
            MappedSparseMatrix_cm A=MappedSparseMatrix_cm(this->height(),this->width(),a_rows.size(),&a_columns[0],&a_rows[0],&(*(data_begin()+(a_temp_begin-a_index_begin)))); //map data to a compressed column matrix

            //compute LU decomposition
            LU_decomp lu_of_A(A);
            if(!lu_of_A.succeeded())
            {
                std::stringstream t;
                t << "Could not perform LU decomp on tab " << k << ".";
                throw RankDeficientException(t.str());
            }

            //now solve for each column of b
            other_storage_iterator cur_begin=b_begin+(b_temp_begin-b_index_begin); //beginning of current tab of b
            other_storage_iterator tab_end=b_begin+(b_temp_end-b_index_begin);  //end of current tab of b

            //loop through every column of b
            for (int j=0; j<b.width(); j++)
            {
                b_vector.setZero(); //reset dense temp vector
                //find upper bound of current column in current tab of b
                other_storage_iterator cur_end=std::upper_bound(cur_begin,tab_end,j,boost::bind(std::less<b_index_type>(),_1, boost::bind(column_other,&b_d,boost::bind(index_bind_other,&b_d,_2))));
                //store nonzero values in dense temp vector
                for (; cur_begin<cur_end; cur_begin++)
                {
                    b_vector(b.row(b.index(*cur_begin)))=to_mdata_type::convert(b.data_val(*cur_begin));

                }
                //get wrapper for corresponding column of lattice c
                c_vector_type c_vector=c.column_vector(j,k);
                //solve and store in lattice c
                if(!lu_of_A.solve(b_vector,&c_vector))
                {
                    std::stringstream t;
                    t << "Solution process on tab " << k << " and column "<< j << "of RHS failed.";
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





/*! @} */

} //libMIA





#endif // SPARSELATTICEBASE_H
