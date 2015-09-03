/*
Taken from Anthony Williams contribution on Boost's Yahoo groups
http://groups.yahoo.com/group/boost/files/tupleit.zip and w
Modified to use std::pair instead of boost::tuple, which sped things up.
APH-
Even so, it was faster to use an implementation of IntroSort and modify it
to swap one array based on the first array
*/


#ifndef ITERATORS_TUPLEIT_HH
#define ITERATORS_TUPLEIT_HH

#include <iterator>
#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include <new>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/utility.hpp"
#include "boost/type_traits.hpp"
#include "boost/optional.hpp" // for aligned_storage
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <memory>

namespace iterators
{
    namespace detail
    {

        template<typename TupleType>
        void preincrementTuple(TupleType& lhs)
        {

            ++(std::get<0>(lhs));
            ++(std::get<1>(lhs));
        }


        template<typename TupleType>
        void predecrementTuple(TupleType& lhs)
        {

            --(std::get<0>(lhs));
            --(std::get<1>(lhs));
        }



        template<typename difference_type,typename TupleType>
        void addToTuple(TupleType& lhs,difference_type diff)
        {

            std::get<0>(lhs)+=diff;
            std::get<1>(lhs)+=diff;
        }



        template<typename difference_type,typename TupleType>
        void subFromTuple(TupleType& lhs,difference_type diff)
        {

            std::get<0>(lhs)-=diff;
            std::get<1>(lhs)-=diff;
        }

        template<typename difference_type,typename TupleType>
        difference_type diffTuples(TupleType const& lhs,TupleType const& rhs);

        template<typename difference_type,typename TupleType>
        struct DiffTupleHelper
        {
            static difference_type doDiff(TupleType const& lhs,TupleType const& rhs)
            {
                return std::get<0>(lhs)-std::get<0>(rhs);

            }
        };


        template<typename difference_type,typename TupleType>
        difference_type diffTuples(TupleType const& lhs,TupleType const& rhs)
        {
            return DiffTupleHelper<difference_type,TupleType>::doDiff(lhs,rhs);
        }



        template<typename SourceTuple,bool const_refs>
        struct MakeTupleTypeWithReferences
        {

            typedef typename boost::mpl::if_c<
                const_refs,
                std::pair<
                     typename boost::add_reference<const typename SourceTuple::first_type>::type,
                     typename boost::add_reference<const typename SourceTuple::second_type>::type
                >,
                std::pair<
                    typename boost::add_reference<typename SourceTuple::first_type>::type,
                    typename boost::add_reference<typename SourceTuple::second_type>::type
                >
            >::type Type;



            template<typename Tuple>
            static Type makeTuple(Tuple& source)
            {
                return Type(std::get<0>(source),std::get<1>(source));
            }
        };



	typedef char Tiny;
	struct Small
	{
	    Tiny dummy[2];
	};
	struct Medium
	{
	    Small dummy[2];
	};
	struct Large
	{
	    Medium dummy[2];
	};
	struct Huge
	{
	    Large dummy[2];
	};

	template<unsigned>
	struct CategoryMap
	{
	    typedef void Type;
	};



//     Tiny categoryCheck(std::output_iterator_tag*);
	Small categoryCheck(std::input_iterator_tag*);
	Medium categoryCheck(std::forward_iterator_tag*);
	Large categoryCheck(std::bidirectional_iterator_tag*);
	Huge categoryCheck(std::random_access_iterator_tag*);


//     template<>
//     struct CategoryMap<sizeof(Tiny)>
//     {
// 	typedef std::output_iterator_tag Type;
//     };

	template<>
	struct CategoryMap<sizeof(Small)>
	{
	    typedef std::input_iterator_tag Type;
	};

	template<>
	struct CategoryMap<sizeof(Medium)>
	{
	    typedef std::forward_iterator_tag Type;
	};

	template<>
	struct CategoryMap<sizeof(Large)>
	{
	    typedef std::bidirectional_iterator_tag Type;
	};
	template<>
	struct CategoryMap<sizeof(Huge)>
	{
	    typedef std::random_access_iterator_tag Type;
	};

	template<typename Cat1,typename Cat2>
	struct CommonCategory
	{
	private:
	    enum
	    {categorySize=sizeof(::iterators::detail::categoryCheck(false?(Cat1*)0:(Cat2*)0))
	    };
	public:
	    typedef typename CategoryMap<categorySize>::Type Type;
	};

	// specializations
	template<typename Cat>
	struct CommonCategory<std::output_iterator_tag,Cat>
	{
	    typedef std::output_iterator_tag Type;
	};
	template<typename Cat>
	struct CommonCategory<Cat,std::output_iterator_tag>
	{
	    typedef std::output_iterator_tag Type;
	};
	template<>
	struct CommonCategory<std::output_iterator_tag,std::output_iterator_tag>
	{
	    typedef std::output_iterator_tag Type;
	};
	template<>
	struct CommonCategory<std::input_iterator_tag,std::output_iterator_tag>
	{
	    // no Type, because error
	};
	template<>
	struct CommonCategory<std::output_iterator_tag,std::input_iterator_tag>
	{
	    // no Type, because error
	};



        template<typename IterTuple,typename SourceTuple>
        void derefAndWrite(IterTuple& iters,SourceTuple const& source)
        {
            *(std::get<0>(iters))=std::get<0>(source);
            *(std::get<1>(iters))=std::get<1>(source);
        }

    }

    // An OutputTuple holds a tuple of references to iterators, and writes to them on assignment
    template<typename IterTuple,bool const_refs>
    struct OutputTuple:
        public detail::MakeTupleTypeWithReferences<IterTuple,const_refs>::Type,
        boost::noncopyable
    {
    private:
        typedef detail::MakeTupleTypeWithReferences<IterTuple,const_refs> BaseTypeBuilder;
        typedef typename BaseTypeBuilder::Type BaseType;
    public:
	OutputTuple(IterTuple& iters):
            BaseType(BaseTypeBuilder::makeTuple(iters))
	{}

	template<typename SomeTuple>
	OutputTuple& operator=(const SomeTuple& other)
	{
            detail::derefAndWrite(static_cast<BaseType&>(*this),other);
	    return *this;
	}
    };

    // An OwningRefTuple holds a tuple of references,
    // which may point to data within the tuple, or external to it

    namespace detail
    {
        struct PreserveReferences
        {};

        template<typename OwnedType>
        struct OwningBase
        {
            std::auto_ptr<OwnedType> tupleBuf;

            OwningBase()
            {}

            template<typename SomeType>
            OwningBase(SomeType &source):
                tupleBuf(new OwnedType(source))
            {}

        };
    }

    template<typename TupleType,bool const_refs>
    struct OwningRefTuple:
        private detail::OwningBase<TupleType>,
        public detail::MakeTupleTypeWithReferences<TupleType,const_refs>::Type
    {
    private:
        typedef detail::MakeTupleTypeWithReferences<TupleType,const_refs> BaseTypeBuilder;
        typedef typename BaseTypeBuilder::Type BaseType;
        typedef detail::OwningBase<TupleType> OwningBaseType;
    public:

        typedef typename BaseType::first_type first_type;
        typedef typename BaseType::second_type second_type;

    private:
        typedef TupleType OwnedTuple;

	OwnedTuple* getTuplePtr()
	{
	    return this->tupleBuf.get();
	}
    public:
	// copy from other types of tuples too
	template<typename SomeTuple>
	OwningRefTuple(const SomeTuple& other):
            OwningBaseType(other),BaseType(BaseTypeBuilder::makeTuple(*getTuplePtr()))
	{
	}
	// copying copies values by default
	OwningRefTuple(const OwningRefTuple& other):
            OwningBaseType(other),BaseType(BaseTypeBuilder::makeTuple(*getTuplePtr()))
	{
	}

	// allow user to specify
	// whether to preserve references
	template<typename SomeTuple>
	OwningRefTuple(SomeTuple& other,detail::PreserveReferences const&):
            BaseType(BaseTypeBuilder::makeTuple(other))
	{
	}

	// assignment assigns to referenced values
	template<typename SomeTuple>
	OwningRefTuple& operator=(const SomeTuple& other)
	{
            BaseType::operator=(other);
	    return *this;
	}
	OwningRefTuple& operator=(const OwningRefTuple& other)
	{
            BaseType::operator=(other);
	    return *this;
	}
    };

    namespace detail
    {
        template<typename IterTuple>
        struct DerefIterTupleHelperKeepRef
        {
            typedef std::pair<typename std::iterator_traits<typename IterTuple::first_type>::reference,
                                        typename std::iterator_traits<typename IterTuple::second_type>::reference> Type;
        };

        template<typename IterTuple>
        struct DerefIterTupleHelperNoRef
        {
            typedef std::pair<typename std::iterator_traits<typename IterTuple::first_type>::value_type,
                                        typename std::iterator_traits<typename IterTuple::second_type>::value_type> Type;
        };



        template<typename IterTuple>
        const typename DerefIterTupleHelperKeepRef<IterTuple>::Type derefIterTupleKeepRef(IterTuple& iters)
        {
            return typename DerefIterTupleHelperKeepRef<IterTuple>::Type(*std::get<0>(iters),*std::get<1>(iters));
        }


        template<typename IterTuple>
        typename DerefIterTupleHelperNoRef<IterTuple>::Type derefIterTupleNoRef(IterTuple& iters)
        {
            return typename DerefIterTupleHelperNoRef<IterTuple>::Type(*std::get<0>(iters),*std::get<1>(iters));
        }

	// Define, construct and destroy the appropriate value_type for
	// the given iterator category
	template<typename Category,typename IterTuple>
	struct ValueForCategory
	{
        private:
            typedef typename IterTuple::first_type FirstIterType;
            typedef typename IterTuple::second_type SecondIterType;
	    typedef typename std::iterator_traits<FirstIterType>::value_type FirstValueType;
	    typedef typename std::iterator_traits<SecondIterType>::value_type SecondValueType;
	    typedef typename std::iterator_traits<FirstIterType>::reference FirstReferenceType;
        static constexpr bool is_const=boost::is_const<
                                typename boost::remove_reference<
                                    FirstReferenceType
                                >::type
                            >::value;
        public:
            typedef std::pair<FirstValueType,SecondValueType> ValueTuple;

	    typedef OwningRefTuple<ValueTuple,is_const> value_type;
	    typedef value_type Type;

	    static void construct(Type* p,IterTuple const& iters)
	    {
		// don't copy values, keep as references
		new (p) Type(derefIterTupleKeepRef(iters),::iterators::detail::PreserveReferences());
	    }

	    static void destruct(Type* p)
	    {
		p->~OwningRefTuple<ValueTuple,is_const>();
	    }
	};


	template<typename IterTuple>
	struct ValueForCategory<std::input_iterator_tag,IterTuple>
	{
        private:
            typedef typename IterTuple::first_type FirstIterType;
            typedef typename IterTuple::second_type SecondIterType;
	    typedef typename std::iterator_traits<FirstIterType>::value_type FirstValueType;
	    typedef typename std::iterator_traits<SecondIterType>::value_type SecondValueType;
	    typedef typename std::iterator_traits<FirstIterType>::reference FirstReferenceType;
        static constexpr bool is_const=boost::is_const<
                                typename boost::remove_reference<
                                    FirstReferenceType
                                >::type
                            >::value;
        public:
            typedef std::pair<FirstValueType,SecondValueType> ValueTuple;

	    typedef OwningRefTuple<ValueTuple,is_const> value_type;
	    typedef value_type Type;

	    static void construct(Type* p,IterTuple const& iters)
	    {
		// copy values
		new (p) Type(derefIterTupleNoRef(iters));
	    }

	    static void destruct(Type* p)
	    {
		p->~OwningRefTuple<ValueTuple,is_const>();
	    }
	};

	template<typename IterTuple>
	struct ValueForCategory<std::output_iterator_tag,IterTuple>
	{
        public:
	    typedef typename IterTuple::first_type FirstIterType;
	    typedef typename std::iterator_traits<FirstIterType>::reference FirstReferenceType;
	    static constexpr bool is_const=boost::is_const<
                                typename boost::remove_reference<
                                    FirstReferenceType
                                >::type
                            >::value;
	    typedef OutputTuple<IterTuple,is_const> value_type;
	    typedef value_type Type;

	    static void construct(Type* p,IterTuple& iters)
	    {
		// copy values
		new (p) Type(iters);
	    }

	    static void destruct(Type* p)
	    {
		p->~OutputTuple<IterTuple,is_const>();
	    }
	};




	template<typename Category,typename IterTuple>
	struct VFCSelector
	{
	    typedef ValueForCategory<Category,IterTuple> Type;
	};

	// Select the iterator_category and value_type for our TupleIt
	template<typename IterTuple>
	struct TupleItHelper
	{
            typedef typename IterTuple::first_type FirstIterType;
            typedef typename IterTuple::second_type SecondIterType;

	    typedef typename std::iterator_traits<FirstIterType>::iterator_category Cat1;
	    typedef typename std::iterator_traits<SecondIterType>::iterator_category Cat2;

	    typedef typename CommonCategory<Cat1,Cat2>::Type iterator_category;
	    typedef typename VFCSelector<iterator_category,IterTuple>::Type ValueTypeDef;
	    typedef typename ValueTypeDef::value_type value_type;
	    typedef typename ValueTypeDef::Type DeRefType;

	    typedef DeRefType& reference;
	    typedef DeRefType* pointer;

	    typedef std::ptrdiff_t difference_type;

	    typedef std::iterator<iterator_category,value_type,difference_type,pointer,reference> IteratorType;

	    static void construct(DeRefType* p,IterTuple& iters)
	    {
		ValueTypeDef::construct(p,iters);
	    }

	    static void destruct(DeRefType* p)
	    {
		ValueTypeDef::destruct(p);
	    }
	};


    }

    // the actual Tuple Iterator itself
    template<typename IterTuple>
    struct TupleIt:
	public detail::TupleItHelper<IterTuple>::IteratorType
    {
    private:
	typedef detail::TupleItHelper<IterTuple> TupleDefs;
    public:
	typedef typename TupleDefs::iterator_category iterator_category;
	typedef typename TupleDefs::value_type value_type;
	typedef typename TupleDefs::difference_type difference_type;
	typedef typename TupleDefs::reference reference;
	typedef typename TupleDefs::pointer pointer;
    private:
	pointer getValuePtr() const
	{
	    return reinterpret_cast<pointer>(dataCache.address());
	}

	void emptyCache() const
	{
	    if(cacheInitialized)
	    {
		TupleDefs::destruct(getValuePtr());
		cacheInitialized=false;
	    }
	}

	void initCache() const
	{
	    emptyCache();
	    TupleDefs::construct(getValuePtr(),iters);
	    cacheInitialized=true;
	}


    public:

	TupleIt(IterTuple iters_):
	    iters(iters_),cacheInitialized(false)
	{}
	template<typename OtherIterTuple>
	TupleIt(const TupleIt<OtherIterTuple>& other):
	    iters(other.iters),cacheInitialized(false)
	{}
	TupleIt(const TupleIt& other):
	    iters(other.iters),cacheInitialized(false)
	{}
	TupleIt():
	    iters(),cacheInitialized(false)
	{}


	~TupleIt()
	{
	    emptyCache();
	}

	void swap(TupleIt& other)
	{
	    using std::swap;

	    swap(iters,other.iters);
	}

        TupleIt& operator=(TupleIt const& other)
        {
            emptyCache();
            iters=other.iters;
            return *this;
        }

	// Input Iterator requirements
	reference operator*() const
	{
	    initCache();
	    return *getValuePtr();
	}

	pointer operator->() const
	{
	    initCache();
	    return getValuePtr();
	}

	friend bool operator==(const TupleIt& lhs,const TupleIt& rhs)
	{
	    return lhs.iters==rhs.iters;
	}

	friend bool operator!=(const TupleIt& lhs,const TupleIt& rhs)
	{
	    return lhs.iters!=rhs.iters;
	}

	// Forward Iterator requirements
	TupleIt& operator++()
	{
            detail::preincrementTuple(iters);
	    return *this;
	}

	TupleIt operator++(int)
	{
	    TupleIt temp(*this);
	    ++*this;
	    return temp;
	}

	// Bidirectional Iterator requirements
	TupleIt& operator--()
	{
            detail::predecrementTuple(iters);
	    return *this;
	}

	TupleIt operator--(int)
	{
	    TupleIt temp(*this);
	    --*this;
	    return temp;
	}

	// Random-Access Iterator requirements
	TupleIt& operator+=(difference_type n)
	{
            detail::addToTuple(iters,n);
	    return *this;
	}

	TupleIt& operator-=(difference_type n)
	{
            detail::subFromTuple(iters,n);
	    return *this;
	}

	friend difference_type operator-(const TupleIt& a,const TupleIt& b)
	{
            return detail::diffTuples<difference_type>(a.iters,b.iters);
	}

	value_type operator[](difference_type n) const
	{
	    return *(*this+n);
	}

    private:
	// everything is mutable so we can modify it without affecting const correctness
	// of client code
	mutable IterTuple iters;
	mutable boost::optional_detail::aligned_storage<typename TupleDefs::DeRefType> dataCache;
	mutable bool cacheInitialized;
    };

    // more random-access iterator requirements
    template<typename IterTuple>
    TupleIt<IterTuple> operator+(std::ptrdiff_t n,TupleIt<IterTuple> temp)
    {
	temp+=n;
	return temp;
    }

    template<typename IterTuple>
    TupleIt<IterTuple> operator+(TupleIt<IterTuple> temp,std::ptrdiff_t n)
    {
	temp+=n;
	return temp;
    }

    template<typename IterTuple>
    TupleIt<IterTuple> operator-(TupleIt<IterTuple> temp,std::ptrdiff_t n)
    {
	temp-=n;
	return temp;
    }

    template<typename IterTuple,typename IterTuple2>
    bool operator<(const TupleIt<IterTuple>& a,const TupleIt<IterTuple2>& b)
    {
	return (b-a)>0;
    }

    template<typename IterTuple,typename IterTuple2>
    bool operator>(const TupleIt<IterTuple>& a,const TupleIt<IterTuple2>& b)
    {
	return b<a;
    }

    template<typename IterTuple,typename IterTuple2>
    bool operator>=(const TupleIt<IterTuple>& a,const TupleIt<IterTuple2>& b)
    {
	return !(b<a);
    }

    template<typename IterTuple,typename IterTuple2>
    bool operator<=(const TupleIt<IterTuple>& a,const TupleIt<IterTuple2>& b)
    {
	return !(b>a);
    }

    // implementation of swap and iter_swap
    template<typename IterTuple>
    void swap(TupleIt<IterTuple>& lhs,TupleIt<IterTuple>& rhs)
    {
	lhs.swap(rhs);
    }

//     template<typename IterTuple,IterTuple2>
//     void iter_swap(const TupleIt<IterTuple>& lhs,const TupleIt<IterTuple2>& rhs)
//     {
// 	lhs.iter_swap(rhs);
//     }

    template<typename Iter1,typename Iter2>
    TupleIt<typename std::pair<Iter1,Iter2> > makeTupleIterator(Iter1 i1,Iter2 i2)
    {
	return TupleIt<typename std::pair<Iter1,Iter2> >(std::make_pair(i1,i2));
    }

//    template<typename Iter1,typename Iter2,typename Iter3>
//    TupleIt<typename boost::tuples::tuple<Iter1,Iter2,Iter3> > makeTupleIterator(Iter1 i1,Iter2 i2,Iter3 i3)
//    {
//	return TupleIt<typename boost::tuples::tuple<Iter1,Iter2,Iter3> >(boost::make_tuple(i1,i2,i3));
//    }
//
//    template<typename Iter1,typename Iter2,typename Iter3,typename Iter4>
//    TupleIt<typename boost::tuples::tuple<Iter1,Iter2,Iter3,Iter4> > makeTupleIterator(Iter1 i1,Iter2 i2,Iter3 i3,Iter4 i4)
//    {
//	return TupleIt<typename boost::tuples::tuple<Iter1,Iter2,Iter3,Iter4> >(boost::make_tuple(i1,i2,i3,i4));
//    }

}

#endif
