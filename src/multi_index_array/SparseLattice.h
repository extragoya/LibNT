#ifndef SPARSE_LATTICE_H_INCLUDED
#define SPARSE_LATTICE_H_INCLUDED

#include <vector>
#include <boost/operators.hpp>
#include <boost/assert.hpp>
#include <boost/array.hpp>

#include "Util.h"
#include "SparseLatticeBase.h"
#include "DenseLattice.h"







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
    typedef typename boost::tuple<typename const_data_iterator<SparseLattice<T> >::type,typename const_index_iterator<SparseLattice<T> >::type> type;

};

template<typename T>
struct full_tuple<SparseLattice<T> >
{
    typedef boost::tuple<T&,typename index_type<SparseLattice<T> >::type &> type;
};

template<typename T>
struct const_full_tuple<SparseLattice<T> >
{
    typedef boost::tuple< const T&,const typename index_type<SparseLattice<T> >::type &> type;
};

template<typename T>
struct storage_iterator<SparseLattice<T> >
{
    typedef typename iterators::TupleIt<typename full_iterator_tuple<SparseLattice<T> >::type > type;
};

template<typename T>
struct const_storage_iterator<SparseLattice<T> >
{
    typedef typename boost::zip_iterator<typename const_full_iterator_tuple<SparseLattice<T> >::type > type;
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

    //!  Constructs a sparse lattice from two prexisting std::vectors.
    /*!
    Swaps the contents of _data and _indices with the lattice's internal vectors, making passed in vectors empty.
    Will throw a LatticeParameterException upon invalid parameters.
    \param[in] _data std::vector<T> of data values.
    \param[in] _indices std::vector<long long> of indice values. Must be linear index of column/row major ordering
    */
    SparseLattice(Data& _data,Indices& _indices,index_type _height,index_type _width,index_type _depth): m_data(), m_indices()
    {
        if (_data.size()!=_indices.size())
            throw LatticeParameterException("Data and Indices vectors must be the same size.");
        std::swap(m_data,_data);
        std::swap(m_indices,_indices);
        this->init(_height,_width,_depth);
        this->sparse_init(false,ColumnMajor);
        //this->sort(ColumnMajor);

    }

    //!  Copy constructor
    SparseLattice(const SparseLattice& other): m_data(other.m_data),m_indices(other.m_indices)
    {
        this->init(other.height(),other.width(),other.depth());
        this->sparse_init(other.is_sorted(),other.sort_order());

    }

    //!  Copy constructor - accepts other types of SparseLattices
    template<class otherDerived>
    SparseLattice(const SparseLatticeBase<otherDerived>& other)
    {


        for(typename internal::const_storage_iterator<otherDerived>::type other_begin=other.derived().begin(); other_begin<other.derived().end();other_begin++){


            push_back(other.data_val(*other_begin),other.index(*other_begin));

        }

        this->init(other.height(),other.width(),other.depth());
        this->sparse_init(other.is_sorted(),other.sort_order());

    }

    //!  Constructs a sparse lattice from two prexisting std::vectors.
    /*!
    Copies the contents of _data and _indices to the lattice's internal vectors.
    Will throw a LatticeParameterException upon invalid parameters.
    \param[in] _data std::vector<T> of data values.
    \param[in] _indices std::vector<long long> of indice values. Must be linear index of column/row major ordering
    */
    SparseLattice(const Data& _data,const Indices& _indices,index_type _height,index_type _width,index_type _depth):m_data(_data), m_indices(_indices)
    {

        if (_data.size()!=_indices.size())
            throw LatticeParameterException("Data and Indices vectors must be the same size.");
        this->init(_height,_width,_depth);
        this->sparse_init(false,ColumnMajor);

    }

    //!  Copy constructor from a dense lattice
    SparseLattice(const DenseLattice<T>& dlat)
    {
        clear();
        this->init(dlat.height(),dlat.width(),dlat.depth());
        this->sparse_init(false,ColumnMajor);
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


        clear();
        std::copy(b.begin(),b.end(),this->begin());
        this->init(b.height(),b.width(),b.depth());
        this->sparse_init(b.is_sorted(),b.sort_order());

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

    index_iterator index_begin() ;
    index_iterator index_end() ;
    data_iterator data_begin();
    data_iterator data_end();
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
    b.sort();


    typedef  typename internal::data_type<otherDerived >::type b_data_type;
    typedef boost::numeric::converter<data_type,b_data_type> to_mdata_type;
    typedef typename internal::storage_iterator<otherDerived >::type other_storage_iterator;
    other_storage_iterator b_begin=b.begin();
    other_storage_iterator b_end=b.end();

    m_data.reserve(m_data.size()+b.size());
    m_indices.reserve(m_data.size()+b.size());
    storage_iterator a_begin=this->begin();
    storage_iterator a_end=this->end();

    while(a_begin<a_end && b_begin<b_end){
        if (idx_less(index(*a_begin),b.index(*b_begin))){
            a_begin++;
        }
        else if  (idx_less(b.index(*b_begin),index(*a_begin))){
            push_back(b.data_val(*b_begin),b.index(*b_begin));
            b_begin++;

        }
        else{
            data_val(*a_begin)=op(data_val(*a_begin),to_mdata_type::convert(b.data_val(*b_begin)));
            a_begin++;
            b_begin++;
        }

    }
    if (a_begin==a_end)
        while (idx_less(b.index(*b_begin),b.index(*b_end))){
            push_back(b.data_val(*b_begin),b.index(*b_begin));
            b_begin++;
        }

    m_data.reserve(m_data.size());
    m_indices.reserve(m_data.size());
    const index_type& (SparseLattice::*index)(const_full_tuple) const = &SparseLattice::index;
    std::inplace_merge(this->begin(),a_end,this->end(),boost::bind(&SparseLattice::idx_less, this,boost::bind(index,this,_1),boost::bind(index,this,_2)));

    return *this;
}





/*! @} */

}// LibMIA


typedef  LibMIA::SparseLattice<double> dSLattice;
typedef  LibMIA::SparseLattice<int> iSLattice;



#endif // SPARSE_DATA_SELECTOR_H_INCLUDED
