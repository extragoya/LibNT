#ifndef MAPPED_SPARSE_LATTICE_H_INCLUDED
#define MAPPED_SPARSE_LATTICE_H_INCLUDED

#include <vector>
#include <boost/operators.hpp>
#include <boost/assert.hpp>
#include <boost/array.hpp>

#include "Util.h"
#include "SparseLatticeBase.h"








namespace LibMIA
{


/** \addtogroup util Utilities
 *  @{
*/

/** \addtogroup lattice_util Lattice Utilities
 *  @{
 */

//traits for MappedSparseLattice
namespace internal {

template<typename T>
struct data_type<MappedSparseLattice<T> >
{
    typedef T type;
};

template<typename T>
struct index_type<MappedSparseLattice<T> >
{
    typedef long long type;
};

template<typename T>
struct Data<MappedSparseLattice<T> >
{
    typedef typename data_type<MappedSparseLattice<T> >::type* type;
};

template<typename T>
struct Indices<MappedSparseLattice<T> >
{
    typedef typename index_type<MappedSparseLattice<T> >::type*  type;
};

template<typename T>
struct index_iterator<MappedSparseLattice<T> >
{
    typedef typename index_type<MappedSparseLattice<T> >::type* type;
};

template<typename T>
struct data_iterator<MappedSparseLattice<T> >
{
    typedef typename data_type<MappedSparseLattice<T> >::type* type;
};

template<typename T>
struct const_index_iterator<MappedSparseLattice<T> >
{
    typedef typename index_type<MappedSparseLattice<T> >::type const* type;
};

template<typename T>
struct const_data_iterator<MappedSparseLattice<T> >
{
    typedef typename data_type<MappedSparseLattice<T> >::type const* type;
};

template<typename T>
struct full_iterator_tuple<MappedSparseLattice<T> >
{
    typedef typename boost::tuple<typename data_iterator<MappedSparseLattice<T> >::type,typename index_iterator<MappedSparseLattice<T> >::type> type;
};

template<typename T>
struct const_full_iterator_tuple<MappedSparseLattice<T> >
{
    typedef typename boost::tuple<typename const_data_iterator<MappedSparseLattice<T> >::type,typename const_index_iterator<MappedSparseLattice<T> >::type> type;

};

template<typename T>
struct full_tuple<MappedSparseLattice<T> >
{
    typedef boost::tuple<T&,typename index_type<MappedSparseLattice<T> >::type &> type;
};

template<typename T>
struct const_full_tuple<MappedSparseLattice<T> >
{
    typedef boost::tuple< const T&,const typename index_type<MappedSparseLattice<T> >::type &> type;
};

template<typename T>
struct storage_iterator<MappedSparseLattice<T> >
{
    typedef typename iterators::TupleIt<typename full_iterator_tuple<MappedSparseLattice<T> >::type > type;
};

template<typename T>
struct const_storage_iterator<MappedSparseLattice<T> >
{
    typedef typename boost::zip_iterator<typename const_full_iterator_tuple<MappedSparseLattice<T> >::type > type;
};



} //namespace internal

/*! @} */
/*! @} */

/** \addtogroup lattice Lattice Classes
 *  @{
*/
//Used to map existing set of data and indices to a sparse lattice. Destruction of data is left to the caller
template<class T>
class MappedSparseLattice: public SparseLatticeBase<MappedSparseLattice<T > >
{

public:


    typedef typename LibMIA::internal::data_type<MappedSparseLattice>::type data_type;
    typedef typename LibMIA::internal::index_type<MappedSparseLattice>::type index_type;
    typedef typename LibMIA::internal::Data<MappedSparseLattice>::type Data;
    typedef typename LibMIA::internal::Indices<MappedSparseLattice>::type Indices;
    typedef typename LibMIA::internal::full_iterator_tuple<MappedSparseLattice>::type full_iterator_tuple;
    typedef typename LibMIA::internal::const_full_iterator_tuple<MappedSparseLattice>::type const_full_iterator_tuple;
    typedef typename LibMIA::internal::full_tuple<MappedSparseLattice>::type full_tuple;
    typedef typename LibMIA::internal::const_full_tuple<MappedSparseLattice>::type const_full_tuple;
    typedef typename LibMIA::internal::storage_iterator<MappedSparseLattice>::type storage_iterator;
    typedef typename LibMIA::internal::const_storage_iterator<MappedSparseLattice>::type const_storage_iterator;
    typedef typename LibMIA::internal::index_iterator<MappedSparseLattice>::type index_iterator;
    typedef typename LibMIA::internal::data_iterator<MappedSparseLattice>::type data_iterator;
    typedef typename LibMIA::internal::const_index_iterator<MappedSparseLattice>::type const_index_iterator;
    typedef typename LibMIA::internal::const_data_iterator<MappedSparseLattice>::type const_data_iterator;


    //!  Constructs a sparse lattice from pre-allocated data.
    /*!
    \param[in] _data Raw pointer of data memory. Must have have been allocated to same size as _size. No checking is done to ensure this requirement is met.
    \param[in] _indices Raw pointer of index memory. Must have have been allocated to same size as _size. No checking is done to ensure this requirement is met. Must be linear index of column/row major ordering
    \param[in] _size Number of nonzeros contained in the sparse lattice. Undefined behaviour if this is different than size of _data and _indices.
    */
    MappedSparseLattice(Data _data,Indices _indices,index_type _size, index_type _height,index_type _width,index_type _depth,bool _sort_order=ColumnMajor):  m_data(_data),m_indices(_indices), m_size(_size)
    {

        this->init(_height,_width,_depth);
        this->sparse_init(false,_sort_order);
        //this->sort(ColumnMajor);

    }





    void clear()
    {

        for (storage_iterator i=this->begin();i<this->end();i++){
            this->data_val(*i)=0;
            this->index(*i)=0;
        }


    }
    index_type size() const
    {
        return m_size;

    }






    index_iterator index_begin() ;
    index_iterator index_end() ;
    data_iterator data_begin();
    data_iterator data_end();
    const_index_iterator index_begin() const;
    const_index_iterator index_end() const ;
    const_data_iterator data_begin() const;
    const_data_iterator data_end() const;
protected:
    Data m_data;
    Indices m_indices;
    index_type m_size;





};






template<typename T>
inline typename MappedSparseLattice<T>::index_iterator
MappedSparseLattice<T>::index_begin()
{


    return &m_indices[0];


}

template<typename T>
inline typename MappedSparseLattice<T>::const_index_iterator
MappedSparseLattice<T>::index_begin() const
{


    return &m_indices[0];


}


template<typename T>
inline typename MappedSparseLattice<T>::index_iterator
MappedSparseLattice<T>::index_end()
{


    return &m_indices[size()];


}

template<typename T>
inline typename MappedSparseLattice<T>::const_index_iterator
MappedSparseLattice<T>::index_end() const
{


    return &m_indices[size()];


}

template<typename T>
inline typename MappedSparseLattice<T>::data_iterator
MappedSparseLattice<T>::data_begin()
{


    return &m_data[0];


}

template<typename T>
inline typename MappedSparseLattice<T>::const_data_iterator
MappedSparseLattice<T>::data_begin() const
{


    return &m_data[0];


}

template<typename T>
inline typename MappedSparseLattice<T>::data_iterator
MappedSparseLattice<T>::data_end()
{


    return &m_data[size()];


}

template<typename T>
inline typename MappedSparseLattice<T>::const_data_iterator
MappedSparseLattice<T>::data_end() const
{


    return &m_data[size()];


}










/*! @} */

}// LibMIA






#endif //MAPPED_SPARSE_LATTICE_H_INCLUDED