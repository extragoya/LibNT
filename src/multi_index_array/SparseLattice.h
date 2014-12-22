// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef SPARSE_LATTICE_H_INCLUDED
#define SPARSE_LATTICE_H_INCLUDED

#include <vector>
#include <boost/operators.hpp>
#include <boost/assert.hpp>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include "LibMIAUtil.h"
#include "SparseLatticeBase.h"
#include "DenseLattice.h"
//#include "tupleit.hh"






namespace LibMIA
{


/** \addtogroup util Utilities
 *  @{
*/

/** \addtogroup lattice_util Lattice Utilities
 *  @{
 */

//traits for SparseLattice
namespace internal {

template<typename T>
struct data_type<SparseLattice<T> >
{
    typedef T type;
};

template<typename T>
struct index_type<SparseLattice<T> >
{
    typedef long long type;
};

template<typename T>
struct Data<SparseLattice<T> >
{
    typedef std::vector<T> type;
};

template<typename T>
struct Indices<SparseLattice<T> >
{
    typedef std::vector<typename index_type<SparseLattice<T> >::type > type;
};

template<typename T>
struct index_iterator<SparseLattice<T> >
{
    typedef typename Indices<SparseLattice<T> >::type::iterator type;
};



template<typename T>
struct data_iterator<SparseLattice<T> >
{
    typedef typename Data<SparseLattice<T> >::type::iterator type;
};

template<typename T>
struct const_index_iterator<SparseLattice<T> >
{
    typedef typename Indices<SparseLattice<T> >::type::const_iterator type;
};

template<typename T>
struct const_data_iterator<SparseLattice<T> >
{
    typedef typename Data<SparseLattice<T> >::type::const_iterator type;
};

template<typename T>
struct full_iterator_tuple<SparseLattice<T> >
{
    typedef typename boost::tuple<typename data_iterator<SparseLattice<T> >::type,typename index_iterator<SparseLattice<T> >::type> type;
};

template<typename T>
struct const_full_iterator_tuple<SparseLattice<T> >
{
	typedef typename boost::tuple<typename const_data_iterator<SparseLattice<T> >::type, typename const_index_iterator<SparseLattice<T> >::type> type;

};

template<typename T>
struct full_tuple<SparseLattice<T> >
{
	typedef boost::tuple<T&, typename index_type<SparseLattice<T> >::type &> type;
};

template<typename T>
struct const_full_tuple<SparseLattice<T> >
{
	typedef boost::tuple< const T&, const typename index_type<SparseLattice<T> >::type &> type;
};

template<typename T>
struct storage_iterator<SparseLattice<T> >
{
    typedef typename boost::zip_iterator<typename full_iterator_tuple<SparseLattice<T> >::type > type;
};

template<typename T>
struct const_storage_iterator<SparseLattice<T> >
{
	typedef typename boost::zip_iterator<typename const_full_iterator_tuple<SparseLattice<T> >::type > type;
};

template<typename T>
struct data_type_ref<SparseLattice<T> >
{
    typedef T & type;
};

template<typename T>
struct const_data_type_ref<SparseLattice<T> >
{
    typedef const T & type;
};



} //namespace internal

/*! @} */
/*! @} */

/** \addtogroup lattice Lattice Classes
 *  @{
*/

//!  Lattice class for sparse data.
/*!
  Supports addition, multiplication, and solution of systems of linear equations.
  Class can be dynamically resized and owns the underlying data and indices.

  \tparam T   the datatype of individual elements.
*/
template<class T>
class SparseLattice: public SparseLatticeBase<SparseLattice<T > >
{

public:


    typedef typename LibMIA::internal::data_type<SparseLattice>::type data_type;
    typedef typename LibMIA::internal::index_type<SparseLattice>::type index_type;
    typedef typename LibMIA::internal::Data<SparseLattice>::type Data;
    typedef typename LibMIA::internal::Indices<SparseLattice>::type Indices;
    typedef typename LibMIA::internal::full_iterator_tuple<SparseLattice>::type full_iterator_tuple;
    typedef typename LibMIA::internal::const_full_iterator_tuple<SparseLattice>::type const_full_iterator_tuple;
    typedef typename LibMIA::internal::full_tuple<SparseLattice>::type full_tuple;
    typedef typename LibMIA::internal::const_full_tuple<SparseLattice>::type const_full_tuple;
    typedef typename LibMIA::internal::storage_iterator<SparseLattice>::type storage_iterator;
    typedef typename LibMIA::internal::const_storage_iterator<SparseLattice>::type const_storage_iterator;
    typedef typename LibMIA::internal::index_iterator<SparseLattice>::type index_iterator;
    typedef typename LibMIA::internal::data_iterator<SparseLattice>::type data_iterator;
    typedef typename LibMIA::internal::const_index_iterator<SparseLattice>::type const_index_iterator;
    typedef typename LibMIA::internal::const_data_iterator<SparseLattice>::type const_data_iterator;


    //!  Constructs an empty sparse lattice.
    SparseLattice(): m_data(), m_indices()
    {
        this->init(0,0,0);
        this->sparse_init(false,ColumnMajor);

    }

    //!  Constructs an empty sparse lattice of specified size.
    SparseLattice(index_type _height,index_type _width,index_type _depth): m_data(), m_indices()
    {
        this->init(_height,_width,_depth);
        this->sparse_init(false,ColumnMajor);

    }




        //!  Constructs a sparse lattice from two prexisting std::vectors.
    /*!
    Copies the contents of _data and _indices to the lattice's internal vectors.
    Will throw a LatticeParameterException upon invalid parameters.
    \param[in] _data std::vector<T> of data values.
    \param[in] _indices std::vector<long long> of indice values. Must be linear index of column/row major ordering
    */
    SparseLattice(const Data& _data,const Indices& _indices,index_type _height,index_type _width,index_type _depth,bool isSorted=false):m_data(_data), m_indices(_indices)
    {

        if (_data.size()!=_indices.size())
            throw LatticeParameterException("Data and Indices vectors must be the same size.");
        this->init(_height,_width,_depth);
        this->sparse_init(isSorted,ColumnMajor);

    }

    //!  Copy constructor
    SparseLattice(const SparseLattice& other): m_data(other.m_data),m_indices(other.m_indices)
    {
        this->init(other.height(),other.width(),other.depth(),other.solveInfo());
        this->sparse_init(other.is_sorted(),other.linIdxSequence());

    }

    //!  Copy constructor - accepts other types of SparseLattices
    template<class otherDerived>
    SparseLattice(const SparseLatticeBase<otherDerived>& other)
    {


        for(typename internal::const_storage_iterator<otherDerived>::type other_begin=other.derived().begin(); other_begin<other.derived().end();other_begin++){


            push_back(other.data_val(*other_begin),other.index_val(*other_begin));

        }

        this->init(other.height(),other.width(),other.depth(),other.solveInfo());
        this->sparse_init(other.is_sorted(),other.linIdxSequence());

    }

    //!  Move constructor
    SparseLattice(SparseLattice&& other): m_data(), m_indices()
    {
        m_data.swap(other.m_data);
        m_indices.swap(other.m_indices);
        this->init(other.height(),other.width(),other.depth(),other.solveInfo());
        this->sparse_init(other.is_sorted(),other.linIdxSequence());

    }

    //!  Move assignment
    SparseLattice& operator=(SparseLattice&& other)
    {
        m_data.swap(other.m_data);
        m_indices.swap(other.m_indices);
        this->init(other.height(),other.width(),other.depth(),other.solveInfo());
        this->sparse_init(other.is_sorted(),other.linIdxSequence());
        return *this;

    }



    //!  Constructs a sparse lattice from two prexisting std::vectors.
    /*!
    Moves the contents of _data and _indices to the lattice's internal vectors.
    Will throw a LatticeParameterException upon invalid parameters.
    \param[in] _data std::vector<T> of data values.
    \param[in] _indices std::vector<long long> of indice values. Must be linear index of column/row major ordering
    */
    SparseLattice( Data&& _data, Indices&& _indices,index_type _height,index_type _width,index_type _depth,bool _isSorted=false):m_data(), m_indices()
    {

        //std::cout << "MOVED sparse lattice" << std::endl;
        m_data.swap(_data);
        m_indices.swap(_indices);
        if (_data.size()!=_indices.size())
            throw LatticeParameterException("Data and Indices vectors must be the same size.");
        this->init(_height,_width,_depth);
        this->sparse_init(_isSorted,ColumnMajor);
        //this->print();

    }

    //!  Copy constructor from a dense lattice
    SparseLattice(const DenseLattice<T>& dlat)
    {

        clear();
        this->init(dlat.height(),dlat.width(),dlat.depth(),dlat.solveInfo());
        this->sparse_init(true,ColumnMajor);
        for (int k=0; k<dlat.depth(); k++)
        {
            for (int j=0; j<dlat.width(); j++)
            {
                for (int i=0; i<dlat.height(); i++)
                {
                    if (T temp=dlat(i,j,k)) //if nonzero
                    {
                        m_indices.push_back(i+j*dlat.height()+k*dlat.height()*dlat.width());
                        m_data.push_back(temp);
                    }


                }
            }
        }


    }

    template<class otherDerived>
    SparseLattice& operator=(const SparseLatticeBase<otherDerived> &b)
    {



        this->resize(b.size());
        std::copy(b.index_begin(),b.index_end(),this->index_begin());
        std::copy(b.data_begin(),b.data_end(),this->data_begin());
        this->init(b.height(),b.width(),b.depth(),b.solveInfo());
        this->sparse_init(b.is_sorted(),b.linIdxSequence());
    }


    SparseLattice& operator=(const SparseLattice &b)
    {



        this->resize(b.size());
        std::copy(b.index_begin(),b.index_end(),this->index_begin());
        std::copy(b.data_begin(),b.data_end(),this->data_begin());
        this->init(b.height(),b.width(),b.depth(),b.solveInfo());
        this->sparse_init(b.is_sorted(),b.linIdxSequence());

        return *this;
    }


    //!  Performs *this+=b.
    /*!
    Scalar operations will be performed using data_type of *this. May result in loss of precision if
    datatype of *this is a subtype of b.
    */
    template<class otherDerived>
    SparseLattice& operator+=(SparseLatticeBase<otherDerived> &b){

        std::plus<T> op;
        return merge(b,op);

    }

    void resize(size_t _size){
        assert(_size<=this->dimesionality());
        this->m_data.resize(_size);
        this->m_indices.resize(_size);
    }

	void reserve(size_t _size){
		this->m_data.reserve(_size);
		this->m_indices.reserve(_size);
	}


    //!  Performs *this-=b.
    /*!
    Scalar operations will be performed using data_type of *this. May result in loss of precision if
    datatype of *this is a subtype of b.
    */
    template<class otherDerived>
    SparseLattice& operator-=(SparseLatticeBase<otherDerived> &b){
        std::minus<T> op;
        return merge(b,op);

    }

    //!  Sets each tab to an identity matrix.
    /*!
        If non-square, it follows <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1MatrixBase.html#a0650b65c6ae6c3d19a138b72a6d68568">this</a> format.
    */
    void eye()
    {

        m_data.clear();
        m_indices.clear();
        index_type eye_dim=std::min(this->height(),this->width());
        m_data.reserve(std::pow(eye_dim,2)*this->depth());
        m_indices.reserve(std::pow(eye_dim,2)*this->depth());
        for(index_type _tab=0;_tab<this->depth();++_tab)
            for(index_type _idx=0;_idx<eye_dim;++_idx)
                this->push_back(1,this->full2lin_index(_idx,_idx,_tab));



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



        this->sort(this->linIdxSequence());
        auto diff=collect_duplicates_function(this->index_begin(), this->index_end(), this->data_begin(), collector);

        resize(diff);
    }

	/*void remove_zeros(){
		decltype(m_data) new_data;
		decltype(m_indices) new_indices;
		new_data.reserve(this->size());
		new_indices.reserve(this->size());
		for (auto it = this->index_begin(); it < this->index_end(); ++it){
			if (std::abs(this->data_at(it))<)
		}
	}*/


    //SparseLattice<T,true> operator*(SparseLattice &b) ;

    void clear()
    {

        m_data.clear();
        m_indices.clear();

    }
    std::size_t size() const
    {
        return m_data.size();

    }

    void push_back(data_type _data,index_type _idx){
        m_data.push_back(_data);
        m_indices.push_back(_idx);
    }

    template<class other_data_type>
    void push_back(other_data_type _data,index_type _idx){
        typedef boost::numeric::converter<data_type,other_data_type> to_mdata_type;
        m_data.push_back(to_mdata_type::convert(_data));
        m_indices.push_back(_idx);
    }

    //void fill(const DenseLattice<T>& dlat);

    const Data& data() const{
        return m_data;
    }

    const Indices& indices() const{
        return m_indices;
    }

    Data& data() {
        return m_data;
    }

    Indices& indices() {
        return m_indices;
    }


	/*data_type  & data_val(storage_iterator it){
		return (*it).get<0>();
	}

	index_type  & index_val(storage_iterator it){
		return (*it).get<1>();
	}

	const data_type  & data_val(storage_iterator it) const {
		return (*it).get<0>();
	}

	const index_type  & index_val(storage_iterator it) const {
		return (*it).get<1>();
	}*/

    index_iterator index_begin() ;
    index_iterator index_end() ;
    data_iterator data_begin();
    data_iterator data_end();

    const_index_iterator index_begin() const ;
    const_index_iterator index_end() const;
    const_data_iterator data_begin()const;
    const_data_iterator data_end()const;

	storage_iterator begin();
	const_storage_iterator begin() const;
	storage_iterator end();
	const_storage_iterator end() const;


protected:
    Data m_data;
    Indices m_indices;

    template <class otherDerived,class Op>
    SparseLattice& merge(SparseLatticeBase<otherDerived> &b, Op op);

    void push_back(const full_tuple a){
        push_back(this->data_val(a),this->index(a));

    }



};


template<typename T>
inline typename SparseLattice<T>::const_storage_iterator
SparseLattice<T>::begin() const
{


	return boost::make_zip_iterator(boost::make_tuple(this->data_begin(), this->index_begin()));


}


template<typename T>
inline typename SparseLattice<T>::const_storage_iterator
SparseLattice<T>::end() const
{


	return boost::make_zip_iterator(boost::make_tuple(this->data_end(), this->index_end()));


}

template<typename T>
inline typename SparseLattice<T>::storage_iterator
SparseLattice<T>::begin()
{


	return boost::make_zip_iterator(boost::make_tuple(this->data_begin(), this->index_begin()));


}


template<typename T>
inline typename SparseLattice<T>::storage_iterator
SparseLattice<T>::end()
{


	return boost::make_zip_iterator(boost::make_tuple(this->data_end(), this->index_end()));


}


template<typename T>
inline typename SparseLattice<T>::const_index_iterator
SparseLattice<T>::index_begin() const
{


    return m_indices.begin();


}


template<typename T>
inline typename SparseLattice<T>::const_index_iterator
SparseLattice<T>::index_end() const
{


    return m_indices.end();


}

template<typename T>
inline typename SparseLattice<T>::const_data_iterator
SparseLattice<T>::data_begin() const
{


    return m_data.begin();


}


template<typename T>
inline typename SparseLattice<T>::const_data_iterator
SparseLattice<T>::data_end() const
{


    return m_data.end();


}




template<typename T>
inline typename SparseLattice<T>::index_iterator
SparseLattice<T>::index_begin()
{


    return m_indices.begin();


}


template<typename T>
inline typename SparseLattice<T>::index_iterator
SparseLattice<T>::index_end()
{


    return m_indices.end();


}

template<typename T>
inline typename SparseLattice<T>::data_iterator
SparseLattice<T>::data_begin()
{


    return m_data.begin();


}


template<typename T>
inline typename SparseLattice<T>::data_iterator
SparseLattice<T>::data_end()
{


    return m_data.end();


}

template <class T>
template <class otherDerived,class Op>
SparseLattice<T>& SparseLattice<T>::merge(SparseLatticeBase<otherDerived> &b, Op op)
{

    otherDerived& b_derived=b.derived();
    this->sort();
    b.sort(this->linIdxSequence());


    //typedef  typename internal::data_type<otherDerived >::type b_data_type;

	SparseLattice<data_type> C(this->height(), this->width(), this->depth());
	internal::outside_merge_sparse_storage_containers(C, *this, b, op); //don't bother doing an inplace_merge, as even std::algorithm uses an out-of-place copy for this
	*this = std::move(C);
	/*typedef typename internal::storage_iterator<otherDerived >::type other_storage_iterator;
    other_storage_iterator b_begin=b.begin();
    other_storage_iterator b_end=b.end();

    m_data.reserve(m_data.size()+b.size());
    m_indices.reserve(m_data.size()+b.size());
    storage_iterator a_begin=this->begin();
    storage_iterator a_end=this->end();

    while(a_begin<a_end && b_begin<b_end){
        if (idx_less(index_val(*a_begin),b.index_val(*b_begin))){
            a_begin++;
        }
		else if (idx_less(b.index_val(*b_begin), index_val(*a_begin))){
			push_back(b.data_val(*b_begin), b.index_val(*b_begin));
            b_begin++;

        }
        else{
            data_val(*a_begin)=op(data_val(*a_begin),to_mdata_type::convert(b.data_val(*b_begin)));
            a_begin++;
            b_begin++;
        }

    }
    if (a_begin==a_end)
	while (idx_less(b.index_val(*b_begin), b.index_val(*b_end))){
		push_back(b.data_val(*b_begin), b.index_val(*b_begin));
            b_begin++;
        }

    m_data.reserve(m_data.size());
    m_indices.reserve(m_data.size());
	const index_type& (SparseLattice::*index_val)(const_full_tuple) const = &SparseLattice::index_val;
    std::inplace_merge(this->begin(),a_end,this->end(),boost::bind(&SparseLattice::idx_less, this,boost::bind(index_val,this,std::placeholders::_1),boost::bind(index_val,this,std::placeholders::_2)));*/

    return *this;
}





/*! @} */



}// LibMIA


typedef  LibMIA::SparseLattice<double> dSLattice;
typedef  LibMIA::SparseLattice<int> iSLattice;



#endif // SPARSE_DATA_SELECTOR_H_INCLUDED
