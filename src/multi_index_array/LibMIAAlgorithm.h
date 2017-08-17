/***************************************************************************
   The code for introsort and natural mergesort is modified from Keith Schwarz's versions.
   Original code can be found online at:
     1. http://www.keithschwarz.com/interesting/code/?dir=introsort
     2. http://www.keithschwarz.com/interesting/code/?dir=natural-mergesort



   1. This is an implementation of introsort, that accepts a user-defined swapper.
   This is necessary for sorting index and data arrays (based on only the index
   arrays), which is a necessary operation in sparse MIAs and lattices. Unlike
   the STL version, this code will only use iter_swap (and not create temp
   buffers), meaning we can safely use distances between iterators to swap entries
   in the data container when we swap indices.

   As well, reliance on STL algorithms was removed, as
   these may or may not use temporary buffers, etc. This required using a non-STL
   implementation of heapsort, as Schawrz's code used the STL version.

   2. natural mergesort, modified along the lines of above

   Keith Schwarz's original file headers are included before each sort implementation

***************************************************************************
 *
 * An implementation of the introsort (introspective sort) algorithm, a
 * hybrid of quicksort, heapsort, and insertion sort that has particularly
 * good runtime behavior.  It is one of the fastest comparison sorting
 * algorithms in use today, and is the usual implementation of the std::sort
 * algorithm provided with the C++ STL.
 *
 * Introsort aims to get the benefits of quicksort (good locality, in-place,
 * fast runtime) without running into any of its degenerate cases.  To do so,
 * the algorithm begins by guessing what the appropriate depth for the
 * quicksort recursion should be, then fires off a quicksort routine.  If
 * quicksort ever makes too many recursive calls, the introsort routine
 * switches over to using heapsort to sort the range.  This means that in
 * the best case, the algorithm runs a standard quicksort with minimal
 * bookkeeping overhead and thus runs extremely quickly.  In the worst case,
 * the algorithm switches to heapsort and avoids the O(n^2) worst-case of
 * quicksort.
 *
 * The algorithm also contains an additional optimization.  Rather than
 * using the O(n lg n) sorts (quicksort and heapsort) to completely sort the
 * input, instead introsort picks some "block size" and then uses the sorts
 * only on subranges larger than the block size.  It then makes a final pass
 * over the input using insertion sort to fix up the range.  Since insertion
 * sort runs extremely quickly (O(n)) when all of the elements in the range
 * are known to be a constant number of positions from their final locations,
 * this step runs rapidly.  It also decreases the overall work necessary by
 * the algorithm, since heapsort and quicksort are expensive on small ranges.
 *
 * This implementation of introsort uses the provided STL implementation of
 * heapsort (make_heap, sort_heap) for simplicity, but has its own versions
 * of the quicksort and insertion sort routines.  It is based on David
 * Musser's original paper on introsort (which can be found at
 * http://www.cs.rpi.edu/~musser/gp/introsort.ps), though it does not use
 * directly any of the code it contains.
 */
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
#ifndef LIBMIAALGORITHM
#define LIBMIAALGORITHM
#define _MIN_GALLOP_ 7
//#include <boost/timer/timer.hpp>
#include <forward_list>
//#include <boost/timer/timer.hpp>
#include "LibMIATimSort.h"
namespace LibMIA
{







namespace internal
{


//little swapper structur for the CO format, just used as part of a benchmark comparison vs the LCO format for sorting
template<typename MainIterator, typename FollowerIterator, size_t co_length>
struct CoSwapper
{
	typedef FollowerIterator FollowerItType;
	typedef MainIterator MainItType;
	const MainIterator _mainBegin;
	const FollowerIterator _followBegin;

	CoSwapper(MainIterator _mainIt, FollowerIterator _followIt) : _mainBegin(_mainIt), _followBegin(_followIt)
	{
	}
	void operator()(MainIterator it1, MainIterator it2) const
	{
		std::iter_swap(_followBegin + (it1 - _mainBegin), _followBegin + (it2 - _mainBegin));
		auto it_raw = it1.m_iter; //now actually get the raw iterators pointing to the packed CO data
		auto it_raw2 = it2.m_iter;



		static_for<0, co_length>::_do([&](int i) //swap values, based on how many numbers are in each CO index
		{

			std::iter_swap(it_raw + i, it_raw2 + i);
		});



	}
	FollowerItType getFollowIt(MainIterator it) const
	{
		return  _followBegin + (it - _mainBegin);

	}

};

//little swapper structur for the CO format as used by MTT, just used as part of a benchmark comparison vs the LCO format for sorting
template<typename MainIterator, typename FollowerIterator, size_t co_length>
struct CoSwapperMTT
{
	typedef FollowerIterator FollowerItType;
	typedef MainIterator MainItType;
	const MainIterator _mainBegin;
	const FollowerIterator _followBegin;

	CoSwapperMTT(MainIterator _mainIt, FollowerIterator _followIt) : _mainBegin(_mainIt), _followBegin(_followIt)
	{
	}
	void operator()(MainIterator it1, MainIterator it2) const
	{
		std::iter_swap(_followBegin + (it1 - _mainBegin), _followBegin + (it2 - _mainBegin));
		auto const it_array = *it1; //now actually get the raw iterators pointing to the packed CO data
		auto const it_array2 = *it2;


		static_for<0, co_length>::_do([&](int i)
		{
			std::iter_swap(it_array[i], it_array2[i]); //swap values, based on how many numbers are in each CO index
		});

	}
	FollowerItType getFollowIt(MainIterator it) const
	{
		return  _followBegin + (it - _mainBegin);

	}

};


//will also swap a second set of data. Two sets of data must correspond with each other.
template<typename MainIterator, typename FollowerIterator>
struct DualSwapper
{
    typedef FollowerIterator FollowerItType;
    typedef MainIterator MainItType;
    const MainIterator _mainBegin;
    const FollowerIterator _followBegin;
    DualSwapper(MainIterator _mainIt,FollowerIterator _followIt): _mainBegin(_mainIt),_followBegin(_followIt)
    {
    }
    void operator()(MainIterator it1, MainIterator it2) const
    {
        std::iter_swap(_followBegin+(it1-_mainBegin),_followBegin+(it2-_mainBegin));
        std::iter_swap(it1,it2);
    }
    FollowerItType getFollowIt(MainIterator it) const
    {
         return  _followBegin+(it-_mainBegin);

    }

};

//will also swap a two more sets of data. Two sets of data must correspond with each other.
template<typename MainIterator, typename FollowerIterator,typename FollowerIterator2>
struct TripleSwapper
{


    const MainIterator _mainBegin;
    const FollowerIterator _followBegin;
    const FollowerIterator2 _followBegin2;
    TripleSwapper(MainIterator _mainIt,FollowerIterator _followIt,FollowerIterator2 _followIt2): _mainBegin(_mainIt),_followBegin(_followIt),_followBegin2(_followIt2)
    {
    }
    void operator()(MainIterator it1, MainIterator it2) const
    {
        std::iter_swap(_followBegin2+(it1-_mainBegin),_followBegin2+(it2-_mainBegin));
        std::iter_swap(_followBegin+(it1-_mainBegin),_followBegin+(it2-_mainBegin));
        std::iter_swap(it1,it2);
    }


};


/**
 * Function: InterpolationSearch(RandomIterator begin, RandomIterator end,
 *                               Element elem);
 * ------------------------------------------------------------------------
 * Performs interpolation search on the sorted range [begin, end).  It is
 * assumed that this range consists of finite integral values and that the
 * input is sorted in ascending order.  Returns whether the element was
 * located.
 */
template <typename RandomIterator, typename Element >
RandomIterator InterpolationSearchUpperBound(RandomIterator begin, RandomIterator end,
                                   const typename std::iterator_traits<RandomIterator>::value_type val,const Element & elemGetter) {


    RandomIterator it;

    typename std::iterator_traits<RandomIterator>::value_type temp;
    typedef typename std::iterator_traits<RandomIterator>::difference_type diffT;
    diffT  step;

    while (begin<end)
    {
        if(elemGetter(*begin)>val)
            return begin;
        it = begin;
        temp=elemGetter(*(end - 1));
        if(temp==val)
            return end;
        //std::cout << "end val is " << elemGetter(*(end - 1)) << " with distance " << end - begin << std::endl;
        if(temp==val+1)
            step=(end-begin)/2;
        else{
            const double interpolation = 1 / (double(temp) - double(val));
            step=diffT(interpolation * (double(end - begin) - 1));
        }
        //std::cout << "step is " << step << std::endl;
        std::advance (it,step);
        //std::cout << "new elem is " << elemGetter(*it) << std::endl;
        if (!(val<elemGetter(*it)))                 // or: if (!comp(val,*it)), for version (2)
        { begin=++it;   }
        else {end=it;
        }

    }

    return begin;



}

template <typename RandomIterator>
RandomIterator InterpolationSearchUpperBound(RandomIterator begin, RandomIterator end,
                                   const typename std::iterator_traits<RandomIterator>::value_type val) {
        RandomIterator it;

    typename std::iterator_traits<RandomIterator>::value_type temp;
    typedef typename std::iterator_traits<RandomIterator>::difference_type diffT;
    diffT  step;

    while (begin<end)
    {
        if(*begin>val)
            return begin;
        it = begin;
        temp=*(end - 1);
        if(temp==val)
            return end;
        //std::cout << "end val is " << elemGetter(*(end - 1)) << " with distance " << end - begin << std::endl;
        if(temp==val+1)
            step=(end-begin)/2;
        else{
            const double interpolation = 1 / (double(temp) - double(val));
            step=diffT(interpolation * (double(end - begin) - 1));
        }
        //std::cout << "step is " << step << std::endl;
        std::advance (it,step);
        //std::cout << "new elem is " << elemGetter(*it) << std::endl;
        if (!(val<*it))                 // or: if (!comp(val,*it)), for version (2)
        { begin=++it;   }
        else {end=it;
        }

    }

    return begin;


}

template <typename RandomIterator>
RandomIterator InterpolationSearchUpperBound(RandomIterator begin, RandomIterator end) {
    return InterpolationSearchUpperBound(begin, end, *begin);


}

/**
 * Function: InterpolationSearch(RandomIterator begin, RandomIterator end,
 *                               Element elem);
 * ------------------------------------------------------------------------
 * Performs interpolation search on the sorted range [begin, end).  It is
 * assumed that this range consists of finite integral values and that the
 * input is sorted in ascending order.  Returns whether the element was
 * located.
 */
template <typename RandomIterator, typename Element >
RandomIterator InterpolationSearchLowerBound(RandomIterator begin, RandomIterator end,
                                   typename std::iterator_traits<RandomIterator>::value_type val,const Element & elemGetter) {


    RandomIterator it;

    typename std::iterator_traits<RandomIterator>::value_type temp;
    typedef typename std::iterator_traits<RandomIterator>::difference_type diffT;
    diffT  step;

    while (begin<end)
    {
        if(elemGetter(*begin)==val)
            return begin;
        it = begin;
        temp=elemGetter(*(end - 1));

        if(temp<=val+1)
            step=(end-begin)/2;
        else{
            const double interpolation = 1 / (double(temp) - double(val));
            step=diffT(interpolation * (double(end - begin) - 1));
        }
        //std::cout << "step is " << step << std::endl;
        std::advance (it,step);
        //std::cout << "new elem is " << elemGetter(*it) << std::endl;
        if (elemGetter(*it)<val)                 // or: if (!comp(val,*it)), for version (2)
        { begin=++it;   }
        else {end=it;
        }

    }


    return begin;



}

/**
 * Function: InterpolationSearch(RandomIterator begin, RandomIterator end,
 *                               Element elem);
 * ------------------------------------------------------------------------
 * Performs interpolation search on the sorted range [begin, end).  It is
 * assumed that this range consists of finite integral values and that the
 * input is sorted in ascending order.  Returns whether the element was
 * located.
 */
template <typename RandomIterator, typename Element >
RandomIterator InterpolationSearchLowerBound(RandomIterator begin, RandomIterator end,
                                   typename std::iterator_traits<RandomIterator>::value_type val) {


    RandomIterator it;

    typename std::iterator_traits<RandomIterator>::value_type temp;
    typedef typename std::iterator_traits<RandomIterator>::difference_type diffT;
    diffT  step;

    while (begin<end)
    {
        if(*begin==val)
            return begin;
        it = begin;
        temp=*(end - 1);

        if(temp<=val+1)
            step=(end-begin)/2;
        else{
            const double interpolation = 1 / (double(temp) - double(val));
            step=diffT(interpolation * (double(end - begin) - 1));
        }

        std::advance (it,step);

        if (*it<val)                 // or: if (!comp(val,*it)), for version (2)
        { begin=++it;   }
        else {end=it;
        }

    }


    return begin;



}


/**
 * Function: Introsort(RandomIterator begin, RandomIterator end);
 * ------------------------------------------------------------------------
 * Sorts the range [begin, end) into ascending order using the introsort
 * algorithm.
 */
template <typename RandomIterator>
void Introsort(RandomIterator begin, RandomIterator end);

/**
 * Function: Introsort(RandomIterator begin, RandomIterator end,
 *                     Comparator comp);
 * -----------------------------------------------------------------------
 * Sorts the range [begin, end) into ascending order (according to comp)
 * using the introsort algorithm.
 */
template <typename RandomIterator, typename Comparator>
void Introsort(RandomIterator begin, RandomIterator end, Comparator comp);

/**
 * Function: Introsort(RandomIterator begin, RandomIterator end,
 *                     Comparator comp);
 * -----------------------------------------------------------------------
 * Sorts the range [begin, end) into ascending order (according to comp)
 * and using swapper to perform the swap using the introsort algorithm.
 */
template <typename RandomIterator, typename Comparator, typename Swapper>
void Introsort(RandomIterator begin, RandomIterator end, Comparator comp, Swapper swapper);

/* * * * * Implementation Below This Point * * * * */
namespace introsort_detail {
  /**
   * Function: Partition(RandomIterator begin, RandomIterator end,
   *                     Comparator comp);
   * Usage: Partition(begin, end, comp);
   * -------------------------------------------------------------
   * Applies the partition algorithm to the range [begin, end),
   * assuming that the pivot element is pointed at by begin.
   * Comparisons are performed using comp.  Returns an iterator
   * to the final position of the pivot element.
   */
  template <typename RandomIterator, typename Comparator, typename Swapper>
  RandomIterator Partition(RandomIterator begin, RandomIterator end,
                           const Comparator & comp,const Swapper & swapper) {
    /* The following algorithm for doing an in-place partition is
     * one of the most efficient partitioning algorithms.  It works
     * by maintaining two pointers, one on the left-hand side of
     * the range and one on the right-hand side, and marching them
     * inward until each one encounters a mismatch.  Once the
     * mismatch is found, the mismatched elements are swapped and
     * the process continues.  When the two endpoints meet, we have
     * found the ultimate location of the pivot.
     */
    RandomIterator lhs = begin + 1;
    RandomIterator rhs = end - 1;
    while (true) {
      /* Keep marching the right-hand side inward until we encounter
       * an element that's too small to be on the left or we hit the
       * left-hand pointer.
       */
      while (lhs < rhs && !comp(*rhs, *begin))
        --rhs;
      /* Keep marching the left-hand side forward until we encounter
       * a the right-hand side or an element that's too big to be
       * on the left-hand side.
       */
      while (lhs < rhs && comp(*lhs,*begin))
        ++lhs;

      /* Now, if the two pointers have hit one another, we've found
       * the crossover point and are done.
       */
      if (lhs == rhs) break;

      /* Otherwise, exchange the elements pointed at by rhs and lhs. */
      swapper(lhs, rhs);
    }
    /* When we have reached this point, the two iterators have crossed
     * and we have the partition point.  However, there is one more edge
     * case to consider.  If the pivot element is the smallest element
     * in the range, then the two pointers will cross over on the first
     * step.  In this case, we don't want to exchange the pivot element
     * and the crossover point.
     */
    if (!comp(*lhs,*begin))
      return begin;

    /* Otherwise, exchange the pivot and crossover, then return the
     * crossover.
     */
    swapper(begin, lhs);
    return lhs;
  }


      /**
   * Function: SiftDown(RandomIterator start, RandomIterator end,RandomIterator current,Comparator comp,Swapper swapper);
   * ---------------------------------------------------------------
   * Sifts element at current
   */
    template <typename RandomIterator, typename Comparator,typename Swapper>
    void SiftDown( RandomIterator begin, ptrdiff_t _size,ptrdiff_t _curIdx,const Comparator & comp,const Swapper &swapper)
    {
        auto root = _curIdx;

        while ( root*2+1 < _size ) {
            auto child = root*2+1;
            if ((child + 1 < _size) && comp(*(begin+child),*(begin+child+1))) {
                child += 1;
            }
            if (comp(*(begin+root), *(begin+child))) {
                swapper(begin+child, begin+root );
                root = child;
            }
            else
                return;
        }
    }

    /**
   * Function: HeapSort(RandomIterator start, RandomIterator end,Comparator comp,Swapper swapper);
   * ---------------------------------------------------------------
   * Main heap sort function
   */
    template <typename RandomIterator, typename Comparator,typename Swapper>
    void HeapSort( RandomIterator begin, RandomIterator end,const Comparator & comp,const Swapper & swapper)
    {
        long long i_start, i_end;
		long long _size = end - begin;
        /* heapify */
        for (i_start = (_size-2)/2; i_start >=0; i_start--) {
            SiftDown(begin,_size,i_start, comp, swapper);
        }

        for (i_end=_size-1; i_end > 0; i_end--) {
            swapper(begin+i_end,begin);
            SiftDown(begin, i_end,0,comp,swapper);
        }
    }




  /**
   * Function: MedianOfThree(RandomIterator one, RandomIterator two,
   *                         RandomIterator three, Comparator comp);
   * ---------------------------------------------------------------
   * Returns the middle element of the three, according to comp.
   */
  template <typename RandomIterator, typename Comparator>
  RandomIterator MedianOfThree(RandomIterator one, RandomIterator two,
                               RandomIterator three, const Comparator& comp) {
    /* Do all three comparisons to determine which is in the middle. */
    const bool comp12 = comp(*one, *two);
    const bool comp13 = comp(*one, *three);
    const bool comp23 = comp(*two, *three);

    /* Based on the relationships between them, return the proper entry. */
    if (comp12 && comp23) return two;               // 1  < 2  < 3
    if (comp12 && !comp23 && comp13) return three;  // 1  < 3 <= 2
    if (!comp12 && comp13) return one;              // 2 <= 1  < 3
    if (!comp12 && !comp13 && comp23) return three; // 2 <  3 <= 1
    if (comp12 && !comp13) return one;              // 3 <= 1  < 2
    return two;                                     // 3 <= 2 <= 1
  }

  /**
   * Function: IntrosortRec(RandomIterator begin, RandomIterator end,
   *                        size_t depth, Comparator comp);
   * ---------------------------------------------------------------------
   * Uses the introsort logic (hybridized quicksort and heapsort) to
   * sort the range [begin, end) into ascending order by comp.
   */
  template <typename RandomIterator, typename Comparator,typename Swapper>
  void IntrosortRec(RandomIterator begin, RandomIterator end,
                    size_t depth, const Comparator & comp,const Swapper & swapper) {
    /* Constant controlling the minimum size of a range to sort.  Increasing
     * this value reduces the amount of recursion performed, but may increase
     * the final runtime by increasing the time it takes insertion sort to
     * fix up the sequence.
     */
    const size_t kBlockSize = 50;


    /* Cache how many elements there are. */
    const size_t numElems = size_t(end - begin);

    /* If there are fewer elements in the range than the block size, we're
     * done.
     */
    if (numElems < kBlockSize) return;

    /* If the depth is zero, sort everything using heapsort, then bail out. */
    if (depth == 0) {
      HeapSort(begin,end,comp,swapper);
      return;
    }

    /* Otherwise, use a median-of-three to pick a (hopefully) good pivot,
     * and partition the input with it.
     */

    RandomIterator pivot = MedianOfThree(begin,                // First elem
                                         begin + numElems / 2, // Middle elem
                                         end - 1, comp);       // Last elem


    /* Swap the pivot in place. */
    swapper(pivot, begin);

    /* Get the partition point and sort both halves. */
    RandomIterator partitionPoint = Partition(begin, end, comp,swapper);
    IntrosortRec(begin, partitionPoint, depth - 1, comp,swapper);
    IntrosortRec(partitionPoint + 1, end, depth - 1, comp,swapper);
  }

  /**
   * Function: IntrosortDepth(RandomIterator begin, RandomIterator end);
   * ---------------------------------------------------------------------
   * Returns the maximum depth to which introsort should be run on a range
   * of the specified size.  This is currently 2 lg (|end - begin|), as
   * suggested in David Musser's paper.
   */
  template <typename RandomIterator>
  size_t IntrosortDepth(RandomIterator begin, RandomIterator end) {
    size_t numElems = size_t(end - begin);

    /* Compute lg(numElems) by shifting the number down until we zero it. */
    size_t lg2 = 0;
    for (; numElems != 0; numElems >>= 1, ++lg2)
      ;

    /* Return twice this value. */
    return lg2 * 2;
  }


}
  /**
   * Function: InsertionSort(RandomIterator begin, RandomIterator end,
   *                         Comparator comp);
   * ----------------------------------------------------------------------
   * Sorts the range [begin, end) into ascending order (according to comp)
   * using insertion sort.
   */
template <typename RandomIterator, typename Comparator, typename Swapper>
inline void  InsertionSort(RandomIterator begin, RandomIterator end,
                     Comparator comp,Swapper swapper) {
    /* Edge case check - if there are no elements or exactly one element,
     * we're done.
     */
    if (begin == end || begin + 1 == end) return;

    /* Starting at the second element and continuing rightward, put each
     * element in its proper position.
     */
    for (RandomIterator itr = begin + 1; itr != end; ++itr) {
      /* Continue swapping down until we hit the beginning or are in the
       * correct position.
       */
      for (RandomIterator test = itr; test != begin && comp(*test, *(test - 1)); --test)
        swapper(test, test - 1);
    }




}

  /**
   * Function: InsertionSortImproved(RandomIterator begin, RandomIterator end,FollowIt followBegin,
   *                         Comparator comp);
   * ----------------------------------------------------------------------
   * Sorts the range [begin, end) into ascending order (according to comp)
   * using insertion sort. Uses improvements found in http://www.drdobbs.com/architecture-and-design/algorithm-improvement-through-performanc/220000504?pgno=1
   */
template <typename RandomIterator, typename Comparator, typename FollowIt>
inline void  InsertionSortImproved(RandomIterator begin, RandomIterator end,FollowIt followBegin,
                     Comparator comp) {


	if (end == begin)
		return;
	auto followIt=followBegin+1;
    for (auto it=begin+1;it<end;++it,++followIt)
    {
        if (comp(*it,*(it-1)))       // no need to do (j > 0) compare for the first iteration
        {
            auto curElement = *it;
            auto followElement = *followIt;
            *it=*(it-1);
            *followIt=*(followIt-1);
            auto it2=it-1;
            auto followIt2=followIt-1;
            for (; it2>begin && curElement < *(it2-1); --it2,--followIt2)
            {
                *it2=*(it2-1);
                *followIt2=*(followIt2-1);

            }
            *it2 = curElement;    // always necessary work/write
            *followIt2 = followElement;    // always necessary work/write
        }
        // Perform no work at all if the first comparison fails - i.e. never assign an element to itself!
    }




}

/**
* Function: InsertionSortImproved(RandomIterator begin, RandomIterator end,FollowIt followBegin,
*                         Comparator comp);
* ----------------------------------------------------------------------
* Sorts the range [begin, end) into ascending order (according to comp)
* using insertion sort. Uses improvements found in http://www.drdobbs.com/architecture-and-design/algorithm-improvement-through-performanc/220000504?pgno=1
*/
template <typename RandomIterator, typename Comparator>
inline void  InsertionSortImproved(RandomIterator begin, RandomIterator end,
	Comparator comp) {


	if (end == begin)
		return;
	
	for (auto it = begin + 1; it<end; ++it)
	{
		if (comp(*it, *(it - 1)))       // no need to do (j > 0) compare for the first iteration
		{
			auto curElement = *it;
			
			*it = *(it - 1);
			
			auto it2 = it - 1;
			
			for (; it2>begin && curElement < *(it2 - 1); --it2)
			{
				*it2 = *(it2 - 1);
				

			}
			*it2 = curElement;    // always necessary work/write
			
		}
		// Perform no work at all if the first comparison fails - i.e. never assign an element to itself!
	}




}


/* Implementation of introsort. */
template <typename RandomIterator, typename Comparator, typename Swapper>
void Introsort(RandomIterator begin, RandomIterator end, Comparator comp, Swapper swapper) {
  /* Give easy access to the utiltiy functions. */
  using namespace introsort_detail;

  /* Fire off a recursive call to introsort using the depth estimate of
   * 2 lg (|end - begin|), as suggested in the original paper.
   */

  IntrosortRec(begin, end, IntrosortDepth(begin, end), comp,swapper);

  /* Use insertion sort to clean everything else up. */
  InsertionSort(begin, end, comp,swapper);
}

/* Non-swapper version calls the swapper version. */
template <typename RandomIterator, typename Comparator>
void Introsort(RandomIterator begin, RandomIterator end, Comparator comp) {
  /* Give easy access to the utiltiy functions. */
  auto swapper=[](RandomIterator one, RandomIterator two){
      std::iter_swap(one,two);
  };
  Introsort(begin,end,comp,swapper);
}

/* Non-comparator version calls the comparator version. */
template <typename RandomIterator>
void Introsort(RandomIterator begin, RandomIterator end) {
  Introsort(begin, end,
            std::less<typename std::iterator_traits<RandomIterator>::value_type>());
}

//tailored made to sort CO format
template <typename RandomIterator, typename FollowIt, size_t co_length>
void IntrosortCO(RandomIterator begin, RandomIterator end, FollowIt followBegin, const std::array<size_t, co_length> & sortOrder) {

	typedef typename std::iterator_traits<RandomIterator>::value_type index_type;

	auto sort_compare = [sortOrder](const index_type& idx1, const index_type& idx2){
		const index_type* it = &idx1;
		const index_type* it2 = &idx2;
		for (int i = co_length - 1; i >= 0; --i){
			auto left = *(it + sortOrder[i]);
			auto right = *(it2 + sortOrder[i]);
			if (left < right){
				return true;
			}
			else if (left > right){
				return false;
			}

		}
		return false;
	};

	typedef CoSwapper<RandomIterator, FollowIt, co_length> Swapper;

	Introsort(begin, end, sort_compare, Swapper(begin, followBegin));


}


//tailored made to sort CO format
template <typename RandomIterator, typename FollowIt, size_t co_length>
void IntrosortCO_MTT(RandomIterator begin, RandomIterator end, FollowIt followBegin, const std::array<size_t, co_length> & sortOrder) {



	typedef typename std::iterator_traits<RandomIterator>::value_type it_value_type;
	auto sort_compare = [sortOrder](const it_value_type& idx_array, const it_value_type& idx_array2){

		bool less_than = true;

		for (int i = co_length - 1; i >= 0; --i){
			auto left = *(idx_array[sortOrder[i]]);
			auto right = *(idx_array2[sortOrder[i]]);
			if (left < right){
				return true;
			}
			else if (left > right){
				return false;
			}

		}
		return false;
	};

	typedef CoSwapperMTT<RandomIterator, FollowIt, co_length> Swapper;

	Introsort(begin, end, sort_compare, Swapper(begin, followBegin));


}

/**************************************************************************
 * File: NaturalMergesort.hh
 * Author: Keith Schwarz (htiek@cs.stanford.edu)
 *
 * An implementation of the natural mergesort algorithm.  Natural mergesort
 * is a bottom-up mergesort algorithm that works by implicitly splitting the
 * input up into a sequence of ascending subranges, then merging adjacent
 * subranges together until the entire input is sorted.  Since at each stage
 * the number of sorted subranges decreases by a factor of two, and there
 * are at most n sorted subranges in the initial input, there can be at most
 * O(lg n) rounds of merging.  Moreover, each merge of two sequences takes
 * at most O(n) time, since the maximum size of any two sequences to merge
 * is at most the size of the input array.  This gives a worst-case upper
 * bound of O(n lg n) runtime.  The merging is done out-of-place for
 * simplicity and thus the algorithm uses O(n) auxiliary storage space.
 *
 * However, this algorithm runs very quickly if the input is already sorted
 * to some initial extent.  In the best case, if the input is already
 * fully-sorted, the algorithm will terminate after running one pass over
 * the input, using only O(n) time.  In this sense, natural mergesort is
 * an adaptive sorting algorithm.
 */

/**
 * Function: NaturalMergesort(ForwardIterator begin, ForwardIterator end);
 * ------------------------------------------------------------------------
 * Sorts the range specified by [begin, end) using the natural mergesort
 * algorithm.  Auxiliary storage space is placed into a temporary vector.
 */
template <typename ForwardIterator>
void NaturalMergesort(ForwardIterator begin, ForwardIterator end);

/**
 * Function: NaturalMergesort(ForwardIterator begin, ForwardIterator end,
 *                            Comparator comp);
 * ------------------------------------------------------------------------
 * Sorts the range specified by [begin, end) using the natural mergesort
 * algorithm.  Auxiliary storage space is placed into a temporary vector,
 * and the sequence is ordered by the strict weak ordering comp.
 */
template <typename ForwardIterator, typename FollowIt,typename Comparator>
void NaturalMergesort(ForwardIterator begin, ForwardIterator end,FollowIt followBegin, std::vector<typename std::iterator_traits<ForwardIterator>::value_type> & mainScratch,
                      std::vector<typename std::iterator_traits<FollowIt>::value_type> & followScratch,
                      Comparator  comp);

/* * * * * Implementation Below This Point * * * * */
namespace naturalmergesort_detail {
  /**
   * Function: SortedRangeEnd(ForwardIterator begin, ForwardIterator end,
   *                          Comparator comp);
   * ---------------------------------------------------------------------
   * Returns an iterator to the end of the longest nondecreasing range
   * starting at begin in the range [begin, end), according to comparison
   * comp.  If the entire sequence is sorted, end is returned.
   */
  template <typename ForwardIterator, typename Comparator>
  ForwardIterator SortedRangeEnd(ForwardIterator begin, ForwardIterator end,
                                 const Comparator & comp) {
    /* Edge case - if the input is empty, we're done. */
    if (begin == end) return end;
//    if(finder(*(end-1))==finder(*begin)) return end;
//    //if(log2((unsigned)(end-begin))< (end-begin)/(finder(*(end-1))-finder(*begin))){
//        InterpolationSearch(begin,end,finder);
//        std::exit(1);
//        return begin;
//
////    }
//    else{
        /* Get an iterator that's one step ahead of begin. */
        ForwardIterator next = begin; ++next;

        /* Keep marching these iterators forward until we find a mismatch or
        * hit the end.  A mismatch occurs when the element after the current
        * is strictly less than the current element.
        */
        for (; next != end && !comp(*next, *begin); ++next, ++begin)
        ;

        /* The above loop stops either when next is the end of the sequence or
        * when next is the crossover point.  In either case, return next.
        */
        return next;
//    }
  }



  /**
   * Function: Merge(ForwardIterator begin, ForwardIterator mid,
   *                 ForwardIterator end, Comparator comp);
   * ---------------------------------------------------------------------
   * Merges the sorted ranges [begin, mid) and [mid, end) together into
   * one sorted range, using comp as the comparator.
   */
  template<typename MainIterator,typename FollowIterator,typename ScratchIterator,typename ScratchIterator2,typename  Comp>
  void Merge(MainIterator begin, MainIterator mid,MainIterator end, FollowIterator followIt,ScratchIterator scratchIt1, ScratchIterator2 scratchIt2,
             const Comp & comp) {



    const size_t kBlockSize = 6;
    if((size_t)(end-begin)<kBlockSize){
        InsertionSort(begin,end,comp,internal::DualSwapper<MainIterator,FollowIterator>(begin,followIt));
        std::copy(begin, end, scratchIt1);
        std::copy(followIt, followIt+(end-begin), scratchIt2);
        return;
    }

    /* Continuously choose the smaller of the two to go in front until some
     * range is consumed.
     */



    MainIterator one = begin, two = mid;
    FollowIterator followOne=followIt, followTwo=followIt+(mid-begin);
    while (one != mid && two != end) {
      if (comp(*one, *two)) { // First sequence has smaller element
        *scratchIt1++=*one;
        *scratchIt2++=*followOne;
        ++one;
        ++followOne;
      } else { // Second sequence has smaller element.
        *scratchIt1++=*two;
        *scratchIt2++=*followTwo;
        ++two;
        ++followTwo;
      }


    }
    //std::cout << " finished initial merge one: " <<one-begin << " two: " << two-begin << std::endl;
    /* Once one of the sequences has been exhausted, one of two cases holds:
     *
     * 1. The first sequence was consumed.  In that case, the rest of the
     *    second sequence is valid and we can just copy the merged sequence
     *    in front of it.
     * 2. The second sequence was consumed.  In that case, we copy the rest
     *    of the first sequence into the merged sequence, then write the
     *    merged sequence back.
     */
    if (two == end){
        std::copy(one, mid, scratchIt1);
        std::copy(followOne, followIt+(mid-begin), scratchIt2);

    }
    else{
        std::copy(two, end, scratchIt1);
        std::copy(followTwo, followIt+(end-begin), scratchIt2);

    }

//    std::copy(scratchSpace1.begin(), sIt1, begin);
//    std::copy(scratchSpace2.begin(), sIt2, swapper.getFollowIt(begin));

  }


template<typename MainIterator,typename FollowIterator,typename ScratchIterator,typename ScratchIterator2,typename  Comp>
bool MergeSortLoop(MainIterator begin, MainIterator end, FollowIterator followBegin, ScratchIterator scratchBegin1, ScratchIterator2 scratchBegin2,const Comp& comp){

    bool haveMerged=false;
    auto scratchIt1=scratchBegin1;
    auto scratchIt2=scratchBegin2;
    auto followIt=followBegin;

    for (auto itr = begin; itr != end; ) {
      /* See how far this range extends. */
      //findtimer.resume();
      auto rangeEnd = SortedRangeEnd(itr, end, comp);
      //findtimer.stop();
      /* If we hit the end of the range, we're done with this iteration. */
      if (rangeEnd == end){
        if(haveMerged){
            std::copy(itr,end,scratchIt1);
            std::copy(followIt,followBegin+(end-begin),scratchIt2);

        }
        break;
      }
      /* See where the end of that range is. */
      //findtimer.resume();
      auto nextRangeEnd = SortedRangeEnd(rangeEnd, end, comp);
      //findtimer.stop();
      /* Merge this range with the range after it. */
      Merge(itr, rangeEnd, nextRangeEnd,followIt, scratchIt1,scratchIt2,comp);

      /* Flag that we did at least one merge so we don't stop early. */
      haveMerged = true;

      /* Advance the iterator to the start of the next sequence, which is
       * directly after the end of the second sorted range.
       */
      itr = nextRangeEnd;
      scratchIt1=scratchBegin1+(itr-begin);
      scratchIt2=scratchBegin2+(itr-begin);
      followIt=followBegin+(itr-begin);
    }


    return haveMerged;


}
}

/* Main implementation of the algorithm. */
template <typename ForwardIterator, typename FollowIt,typename Comparator>
void NaturalMergesort(ForwardIterator begin, ForwardIterator end,FollowIt followBegin, std::vector<typename std::iterator_traits<ForwardIterator>::value_type> & mainScratch,
                      std::vector<typename std::iterator_traits<FollowIt>::value_type> & followScratch,
                      Comparator  comp) {
  /* Make utility functions implicitly available. */
  using namespace naturalmergesort_detail;

  const size_t kBlockSize = 24;
  /* As an edge case, if the input range is empty, we're trivially done. */
  if (end-begin<2) return;

  if((size_t)(end-begin)<kBlockSize){
    InsertionSort(begin,end,comp,internal::DualSwapper<ForwardIterator,FollowIt>(begin,followBegin));
    return;
  }
  /* Determine the type of the element being iterated over. */

    mainScratch.resize(end-begin);
    followScratch.resize(end-begin);



    /* Create a vector of Ts that will hold the merged sequence. */

  /* Track whether the current iteration of the algorithm has made any
   * changes.  If it didn't, then we can finish early.
   */
  bool haveMerged;


  /* Continuously loop, merging together ranges until the input is fully
   * sorted.
   */
   //boost::timer::cpu_timer findtimer;
  size_t ctr=0;
  do {
    /* We have not yet merged anything, so clear this flag. */

    //std::cout << "Iteration " << ctr++ << std::endl;
    /* Scan forward in the loop, looking for adjacent pairs of sorted ranges
     * and merging them.
     */
     //we alternate whether we use the original container or the scrathspace to hold the merged lists
    if(ctr%2==0){
        haveMerged=MergeSortLoop(begin, end, followBegin, mainScratch.begin(), followScratch.begin(),comp);
    }
    else
        haveMerged=MergeSortLoop(mainScratch.begin(), mainScratch.end(), followScratch.begin(), begin, followBegin,comp);


    if(!haveMerged && ctr%2){ //if we are done, and the scratchspace holds the final merged lists, we need a final copy back
        std::copy(mainScratch.begin(), mainScratch.end(),begin);
        std::copy(followScratch.begin(), followScratch.end(),followBegin);
    }
    ctr++;
  } while (haveMerged);
  //std::cout << "finding time in sort " << boost::timer::format(findtimer.elapsed()) << std::endl;

}

template <typename Iter,typename ref_t>
ptrdiff_t gallopLeft(ref_t key, Iter const base, ptrdiff_t const len) {
        assert( len > 0 );
        if (key<=*base)
            return 0;
        typedef ptrdiff_t diff_t;
        diff_t lastOfs = 0;
        diff_t ofs = 1;


        diff_t const maxOfs = len ;
        while(ofs < maxOfs && key > *(base +  ofs)) {
            lastOfs = ofs;
            ofs     = (ofs << 1) + 1;

            if(ofs <= 0) { // int overflow
                ofs = maxOfs;
            }
        }
        if(ofs > maxOfs) {
            ofs = maxOfs;
        }



        assert( -1 <= lastOfs && lastOfs < ofs && ofs <= len );

        return std::lower_bound(base+(lastOfs+1), base+ofs, key) - base;
    }

    template <typename Iter, typename ref_t>
    ptrdiff_t gallopRight(ref_t key, Iter const base, ptrdiff_t const len) {
        assert( len > 0);
        if (key<*base)
            return 0;
        typedef ptrdiff_t diff_t;
        diff_t ofs = 1;
        diff_t lastOfs = 0;


        diff_t const maxOfs = len;
        while(ofs < maxOfs && key>= *(base + ofs)) {
            lastOfs = ofs;
            ofs     = (ofs << 1) + 1;

            if(ofs <= 0) { // int overflow
                ofs = maxOfs;
            }
        }
        if(ofs > maxOfs) {
            ofs = maxOfs;
        }



        assert( 0 <= lastOfs && lastOfs < ofs && ofs <= len );

        return std::upper_bound(base+(lastOfs+1), base+ofs, key) - base;
    }

    template <typename Iter, typename ref_t,typename Comp>
    ptrdiff_t gallopRight(ref_t key, Iter const base, ptrdiff_t const len,Comp comp) {
        assert( len > 0);
        if (comp(key,*base))
            return 0;
        typedef ptrdiff_t diff_t;
        diff_t ofs = 1;
        diff_t lastOfs = 0;


        diff_t const maxOfs = len;
        while(ofs < maxOfs && !comp(key,*(base + ofs))) {
            lastOfs = ofs;
            ofs     = (ofs << 1) + 1;

            if(ofs <= 0) { // int overflow
                ofs = maxOfs;
            }
        }
        if(ofs > maxOfs) {
            ofs = maxOfs;
        }



        assert( 0 <= lastOfs && lastOfs < ofs && ofs <= len );

        return std::upper_bound(base+(lastOfs+1), base+ofs, key,comp) - base;
    }


template <typename Iter,typename ref_t>
ptrdiff_t gallopLeftBack(ref_t key, Iter const base, ptrdiff_t const len) {
        assert( len > 0 );
        if(key>*(base+len-1))
            return len;
        typedef ptrdiff_t diff_t;
        diff_t lastOfs = 0;
        diff_t ofs = 1;
        diff_t hint=len-1;

        diff_t const maxOfs = len;
        while(ofs < maxOfs && key<= *(base + (hint - ofs))) {
            lastOfs = ofs;
            ofs     = (ofs << 1) + 1;

            if(ofs <= 0) {
                ofs = maxOfs;
            }
        }
        if(ofs > maxOfs) {
            ofs = maxOfs;
        }

        diff_t const tmp = lastOfs;
        lastOfs          = hint - ofs;
        ofs              = hint - tmp;

        assert( -1 <= lastOfs && lastOfs < ofs && ofs <= len );

        return std::lower_bound(base+(lastOfs+1), base+ofs, key) - base;
}

template<typename It, typename FollowIt>
void GallopMerge(It const base1, FollowIt const baseFollow1,ptrdiff_t len1, It const base2, FollowIt const baseFollow2,ptrdiff_t len2,It dest,FollowIt destFollow,bool do_copy=true) {
        assert( len1 > 0 && len2 > 0 );

//        std::cout << "First" << std::endl;
//        for(auto it=base1;it<base1+len1;++it)
//            std::cout << *it << " " ;
//        std::cout << std::endl;
//
//        std::cout << "Second" << std::endl;
//        for(auto it=base2;it<base2+len2;++it)
//            std::cout << *it << " " ;
//        std::cout << std::endl;

        typedef typename std::iterator_traits<It>::value_type index_type;
        const size_t kBlockSize = 500;
        if(len1+len2<kBlockSize){
            std::copy(base1, base1+len1, dest);
            std::copy(baseFollow1, baseFollow1+len1, destFollow);

            if(do_copy){
                std::copy(base2, base2+len2, dest+len1);
                std::copy(baseFollow2, baseFollow2+len2, destFollow+len1);
            }
           Introsort(dest, dest+len1+len2, std::less<index_type>(),internal::DualSwapper<It,FollowIt>(dest,destFollow));
            return;

        }

        typedef It iter_t;
        typedef FollowIt iter_t2;
        typedef ptrdiff_t diff_t;
        iter_t cursor1 = base1;
        iter_t cursor2     = base2;


        iter_t2 cursorFollow1 = baseFollow1;
        iter_t2 cursorFollow2     = baseFollow2;

        diff_t k = gallopRight(*base2, base1, len1);
        //std::cout << "First gallop right " << k << std::endl;
        assert( k >= 0 );
        if(k>0){
            dest=std::copy(base1,base1+k,dest);
            destFollow=std::copy(baseFollow1,baseFollow1+k,destFollow);
        }

        cursor1 += k;
        cursorFollow1+=k;
        len1  -= k;
        if(len1 == 0) {
            if (do_copy){
                std::copy(cursor2, cursor2 + len2, dest);
                std::copy(cursorFollow2, cursorFollow2 + len2, destFollow);
            }
            return;
        }


        k = gallopLeftBack(*(cursor1 + (len1 - 1)), base2, len2);
        //std::cout << "First gallop left " << k << " " << len2 << std::endl;
        if(k<len2){
            if(do_copy){
                std::copy_backward(base2+k,base2+len2,dest+len1+len2);
                std::copy_backward(baseFollow2+k,baseFollow2+len2,destFollow+len1+len2);
            }
            len2=k;
        }
        assert(len2>=1);
        *(dest++) = *(cursor2++);
        *(destFollow++) = *(cursorFollow2++);
        if(--len2 == 0) {
            std::copy(cursor1, cursor1 + len1, dest);
            std::copy(cursorFollow1, cursorFollow1 + len1, destFollow);
            return;
        }

        if(len1 == 1) {
            std::copy(cursor2, cursor2 + len2, dest);
            std::copy(cursorFollow2, cursorFollow2 + len2, destFollow);
            *(dest + len2) = *cursor1;
            *(destFollow + len2) = *cursorFollow1;
            return;
        }




        assert( len2 > 0 );
        assert( len1 > 1 );



        //std::cout << "Did setup" << std::endl;


        int minGallop(_MIN_GALLOP_);

        // outer:
        while(true) {
            int count1 = 0;
            int count2 = 0;
            //std::cout << "loop" << std::endl;
            bool break_outer = false;
            do {
                assert( len1 > 1 && len2 > 0 );

                if(*cursor2 < *cursor1) {
                    *(dest++) = *(cursor2++);
                    *(destFollow++) = *(cursorFollow2++);
                    ++count2;
                    count1 = 0;
                    if(--len2 == 0) {
                        break_outer = true;
                        break;
                    }
                }
                else {
                    *(dest++) = *(cursor1++);
                    *(destFollow++) = *(cursorFollow1++);
                    ++count1;
                    count2 = 0;
                    if(--len1 == 1) {
                        break_outer = true;
                        break;
                    }
                }
            } while( (count1 | count2) < minGallop );
            if(break_outer) {
                break;
            }

            do {
                assert( len1 > 1 && len2 > 0 );

                count1 = gallopRight(*cursor2, cursor1, len1);
                assert(count1<len1);
                if(count1 != 0) {
                    std::copy_backward(cursor1, cursor1 + count1, dest + count1);
                    std::copy_backward(cursorFollow1, cursorFollow1 + count1, destFollow + count1);
                    dest    += count1;
                    cursor1 += count1;
                    destFollow    += count1;
                    cursorFollow1 += count1;
                    len1    -= count1;

                    if(len1 <= 1) {
                        break_outer = true;
                        break;
                    }
                }
                *(dest++) = *(cursor2++);
                *(destFollow++) = *(cursorFollow2++);
                if(--len2 == 0) {
                    break_outer = true;
                    break;
                }

                count2 = gallopLeft(*cursor1, cursor2, len2);
                assert(count2<=len2);
                if(count2 != 0) {
                    std::copy(cursor2, cursor2 + count2, dest);
                    std::copy(cursorFollow2, cursorFollow2 + count2, destFollow);
                    dest    += count2;
                    cursor2 += count2;
                    destFollow    += count2;
                    cursorFollow2 += count2;
                    len2    -= count2;
                    if(len2 == 0) {
                        break_outer = true;
                        break;
                    }
                }
                *(dest++) = *(cursor1++);
                *(destFollow++) = *(cursorFollow1++);
                if(--len1 == 1) {
                    break_outer = true;
                    break;
                }

                --minGallop;
            } while( (count1 >= _MIN_GALLOP_) | (count2 >= _MIN_GALLOP_) );
            if(break_outer) {
                break;
            }

            if(minGallop < 0) {
                minGallop = 0;
            }
            minGallop += 2;
        } // end of "outer" loop

        //std::cout << "Finished loop" << std::endl;

        if(len1 == 1) {
            assert( len2 > 0 );
            //std::cout << "len2" << len2 << std::endl;
            std::copy(cursor2, cursor2 + len2, dest);
            std::copy(cursorFollow2, cursorFollow2 + len2, destFollow);
            *(dest + len2) = *cursor1;
            *(destFollow + len2) = *cursorFollow1;

        }
        else {
            assert( len1 != 0 && "Comparision function violates its general contract");
            assert( len2 == 0 );
            assert( len1 > 1 );
            //std::cout << "len1" << len1 << std::endl;
            std::copy(cursor1, cursor1 + len1, dest);
            std::copy(cursorFollow1, cursorFollow1 + len1, destFollow);
        }
        //std::cout << "Finished merge" << std::endl;
}

template<typename ForwardIt, typename FollowIt>
void BlindMerge(ForwardIt begin1, ForwardIt end1, FollowIt followBegin1, ForwardIt begin2, ForwardIt end2,FollowIt followBegin2,ForwardIt output, FollowIt followOutput,int array_location_2,int array_location_output){

//    std::cout << "First" << std::endl;
//    for(auto it=begin1;it<end1;++it)
//        std::cout << *it << std::endl;
//    std::cout << std::endl;
//     std::cout << "Second" << std::endl;
//    for(auto it=begin2;it<end2;++it)
//        std::cout << *it << std::endl;
//    std::cout << std::endl;
    //auto old_output=output;
    //auto total=end1-begin1+end2-begin2;
    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;
    const size_t kBlockSize = 500;
    if(end1-begin1+end2-begin2<kBlockSize){
        auto output_end=std::copy(begin1, end1, output);
        auto it2=std::copy(followBegin1, followBegin1+(end1-begin1), followOutput);

        if(array_location_2 !=array_location_output){
            output_end=std::copy(begin2, end2, output_end);
            std::copy(followBegin2, followBegin2+(end2-begin2), it2);
        }
        Introsort(output, output_end, std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(output,followOutput));
        return;

    }

    if(array_location_2 ==array_location_output){
        while(end2>=begin2 && *(end2-1)>*(end1-1)){
            end2--;
        }
    }


    int inner_loop_iterations=std::min(end1-begin1,end2-begin2);;
    while(inner_loop_iterations){

        int idx=0;
        for(;idx< ROUND_DOWN(inner_loop_iterations, 4);idx+=4){
            if (*begin1<*begin2) { // First sequence has smaller element
                *output++=*begin1++;
                *followOutput++=*followBegin1++;
            } else { // Second sequence has smaller element.
                *output++=*begin2++;
                *followOutput++=*followBegin2++;
            }

            if (*begin1<*begin2) { // First sequence has smaller element
                *output++=*begin1++;
                *followOutput++=*followBegin1++;
            } else { // Second sequence has smaller element.
                *output++=*begin2++;
                *followOutput++=*followBegin2++;
            }

            if (*begin1<*begin2) { // First sequence has smaller element
                *output++=*begin1++;
                *followOutput++=*followBegin1++;
            } else { // Second sequence has smaller element.
                *output++=*begin2++;
                *followOutput++=*followBegin2++;
            }

            if (*begin1<*begin2) { // First sequence has smaller element
                *output++=*begin1++;
                *followOutput++=*followBegin1++;
            } else { // Second sequence has smaller element.
                *output++=*begin2++;
                *followOutput++=*followBegin2++;
            }
        }
        for(;idx< inner_loop_iterations;++idx){
            if (*begin1<*begin2) { // First sequence has smaller element
                *output++=*begin1++;
                *followOutput++=*followBegin1++;
            } else { // Second sequence has smaller element.
                *output++=*begin2++;
                *followOutput++=*followBegin2++;
            }
        }
        inner_loop_iterations=std::min(end1-begin1,end2-begin2);
    }


    if (begin2 == end2){
        std::copy(begin1, end1, output);
        std::copy(followBegin1, followBegin1+(end1-begin1), followOutput);

    }
    else{
        if(array_location_2!=array_location_output){
            std::copy(begin2, end2, output);
            std::copy(followBegin2, followBegin2+(end2-begin2), followOutput);
        }

    }

//    for(auto it=old_output;it<old_output+total;++it)
//        std::cout << *it << " ";
//    std::cout << std::endl;

}

template<typename index_type>
int getNextValidRun(int start, const std::vector<index_type>& array_locations){

    int next_loc;
    for(next_loc=start+1;next_loc<array_locations.size();++next_loc){
        if(array_locations[next_loc]!=0)
            break;
    }
    return next_loc;

}

template<typename ForwardIt, typename FollowIt,typename index_type>
int PingPongMerge(ForwardIt begin1, FollowIt followBegin1, ForwardIt begin2, FollowIt followBegin2,std::vector<index_type> & run_starts, std::vector<index_type> & run_sizes){


    std::vector<index_type> array_locations(run_starts.size());



//    std::cout << "Run starts" << std::endl;
//    for(auto &x: run_starts)
//        std::cout << x << " ";
//    std::cout << std::endl;

    for (auto& x: array_locations){
        x=1;
    }

    int run_idx=0;
    int run_idx_next=1;
    int next_to_begin=1;
    int run_total=run_starts.size();
    auto elem_list_size=run_starts.size();
    ForwardIt runStart2, runEnd2;
    FollowIt runFollowStart2;
    //boost::timer::cpu_timer blind_merge;
    //blind_merge.stop();
    //std::cout << "run totals " << run_total << std::endl;
    while(elem_list_size>1){
        //std::cout << "SANE" << std::endl;
        if(run_idx_next>=run_total||run_sizes[run_idx]+run_sizes[run_idx_next]>run_sizes[0]+run_sizes[next_to_begin]){
            run_idx=0;
            run_idx_next=next_to_begin;
            //std::cout << "We need to start loop over " << run_idx << " " << run_idx_next << std::endl;

        }
        else{
            //assert(run_idx<run_total);
            //assert(run_idx_next<run_total);

            if(array_locations[run_idx_next]==1){
                runStart2=begin1+run_starts[run_idx_next];
                runFollowStart2=followBegin1+run_starts[run_idx_next];
            }
            else{
                runStart2=begin2+run_starts[run_idx_next];
                runFollowStart2=followBegin2+run_starts[run_idx_next];
            }
            runEnd2=runStart2+run_sizes[run_idx_next];
            //std::cout << "About to merge " << run_idx << " and " << run_idx_next << std::endl;
            //std::cout << "start is " << run_starts[run_idx] << " and " << run_starts[run_idx_next] << std::endl;
            //std::cout << "size is " << run_sizes[run_idx] << " and " << run_sizes[run_idx_next] << std::endl;

            if(array_locations[run_idx]==1){

                //blind_merge.resume();
                //GallopMerge(begin1+run_starts[run_idx], followBegin1+run_starts[run_idx],run_sizes[run_idx], runStart2, runFollowStart2,run_sizes[run_idx_next],begin2+run_starts[run_idx],followBegin2+run_starts[run_idx],array_locations[run_idx_next]==1);
                BlindMerge(begin1+run_starts[run_idx],begin1+run_starts[run_idx]+run_sizes[run_idx],followBegin1+run_starts[run_idx],runStart2,runEnd2,runFollowStart2,begin2+run_starts[run_idx],followBegin2+run_starts[run_idx],array_locations[run_idx_next],2);
                //blind_merge.stop();
                array_locations[run_idx]=2;

            }
            else{

               // blind_merge.resume();
                //GallopMerge(begin2+run_starts[run_idx], followBegin2+run_starts[run_idx],run_sizes[run_idx], runStart2, runFollowStart2,run_sizes[run_idx_next],begin1+run_starts[run_idx],followBegin1+run_starts[run_idx],array_locations[run_idx_next]==2);
                BlindMerge(begin2+run_starts[run_idx],begin2+run_starts[run_idx]+run_sizes[run_idx],followBegin2+run_starts[run_idx],runStart2,runEnd2,runFollowStart2,begin1+run_starts[run_idx],followBegin1+run_starts[run_idx],array_locations[run_idx_next],1);
                //blind_merge.stop();
                array_locations[run_idx]=1;

            }
            array_locations[run_idx_next]=0;
            //std::cout <<"Done merge " <<std::endl;

            run_sizes[run_idx]+=run_sizes[run_idx_next];

            run_idx=getNextValidRun(run_idx,array_locations);


            //std::cout << "sanity " << run_idx_next << " and " << next_to_begin << std::endl;
            if (run_idx_next==next_to_begin){
                next_to_begin=run_idx;
                //std::cout << "Next to begin updated " << next_to_begin << std::endl;
            }
            //std::cout << "Got valid run " << run_idx << std::endl;

            run_idx_next=getNextValidRun(run_idx,array_locations);

            //std::cout << "Got next valid run " << run_idx_next << std::endl;

            elem_list_size--;
        }



    }

   // std::cout << "Blind Merge " << boost::timer::format(blind_merge.elapsed()) << std::endl;
    //std::cout << "Finsihed" << elem_list_size << std::endl;
    return array_locations[0];

}

template<typename index_type>
void SortRunsBySize(std::vector<index_type>& run_sizes,std::vector<index_type>& head_offsets,std::vector<index_type>& run_designations){

    std::vector<index_type> run_refs(run_sizes.size());
    std::vector<index_type> run_refs2(run_sizes.size());
    typedef decltype(run_refs.begin()) it_type;


    for(auto _idx=0;_idx<run_refs.size();++_idx){
        run_refs[_idx]=_idx;
    }
//    std::cout << "Run sizes" << std::endl;
//    for(auto &x: run_sizes)
//        std::cout << x << " ";
//    std::cout << std::endl;


    //sort runs by size
    internal::Introsort(run_sizes.begin(),run_sizes.end(),std::less<index_type>(),internal::TripleSwapper<it_type,it_type,it_type>(run_sizes.begin(),head_offsets.begin(),run_refs.begin()));
    //reverse the run specifications
    for(auto _idx=0;_idx<run_refs.size();++_idx){
        run_refs2[run_refs[_idx]]=_idx+1; //run designations are 1-indexed
    }
//    std::cout << "Run refs" << std::endl;
//    for(auto &x: run_refs2)
//        std::cout << x << " ";
//    std::cout << std::endl;
    //update the run designations to the sorted layout
    for(auto it=run_designations.begin();it<run_designations.end();++it){
        if(*it>0)
            *it=run_refs2[*it-1];
        else
            *it=-1*(run_refs2[std::abs(*it)-1]);
    }

}

template <typename ForwardIt,typename FollowIt>
void PerformPatienceRun(ForwardIt begin, ForwardIt end,FollowIt followBegin, std::vector<typename std::iterator_traits<ForwardIt>::value_type> & mainScratch,
                      std::vector<typename std::iterator_traits<FollowIt>::value_type> & followScratch) {


    //boost::timer::cpu_timer total, setup;
    //first let's record the number of runs and their size
    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;
    std::vector<index_type> tails;
    std::vector<index_type> heads;
    std::vector<index_type> run_sizes;
    std::vector<index_type> head_offsets;
    std::vector<index_type> temp_starts;
    std::vector<index_type> run_starts;
    std::vector<index_type> copy_space(end-begin);
    auto cur_it=begin;
    heads.reserve(std::floor(std::sqrt(end-begin)));
    tails.reserve(std::floor(std::sqrt(end-begin)));
    run_sizes.reserve(std::floor(std::sqrt(end-begin)));
    head_offsets.reserve(std::floor(std::sqrt(end-begin)));
    head_offsets.push_back(0);
    heads.push_back(*cur_it);
    index_type cur_min=*cur_it;
    index_type cur_max=*cur_it;
    tails.push_back(*cur_it++);

    run_sizes.push_back(1);
    auto search_region=tails.size();
    copy_space[0]=1;
    decltype(tails.begin()) find_it;
    //boost::timer::cpu_timer creating_runs,finding_timer;
    //std::cout << "Setup " << boost::timer::format(setup.elapsed()) << std::endl;
    //std::cout << "Current tail " << *(tails.begin()) << std::endl;
    while(cur_it<end){
         if(*cur_it>cur_min){
            //finding_timer.resume();
            find_it=std::upper_bound(tails.end()-search_region,tails.end(),*cur_it,std::greater<index_type>());
            //auto k=gallopRight(*cur_it, tails.end()-search_region, search_region,std::greater<index_type>());
            //find_it=tails.end()-search_region+k;
            //finding_timer.stop();
            //std::cout << "for val " << *cur_it << " found tail " << *find_it << std::endl;
            auto run_number=find_it-tails.begin();
            copy_space[cur_it-begin]=run_number+1; //record the run number where the current index should reside in
            *find_it=*cur_it++;
            run_sizes[run_number]++;
            if(find_it!=tails.begin()){
                while(cur_it< end && *cur_it>*find_it && *cur_it<*(find_it-1)){
                    copy_space[cur_it-begin]=run_number+1;
                    *find_it=*cur_it++;
                    run_sizes[run_number]++;

                }
            }
            else{
                while(cur_it<end && *cur_it>*find_it){
                    copy_space[cur_it-begin]=run_number+1;
                    *find_it=*cur_it++;
                    run_sizes[run_number]++;

                }
            }
            if(find_it==tails.end()-1)
                cur_min=*(cur_it-1);
        }
//        else if(*cur_it<cur_max){
//            //finding_timer.resume();
//            //find_it=std::lower_bound(heads.end()-search_region,heads.end(),*cur_it);
//            auto k=gallopRight(*cur_it, heads.end()-search_region, search_region);
//            find_it=heads.end()-search_region+k;
//            //finding_timer.stop();
//            //std::cout << "for val " << *cur_it << " found head " << *find_it << std::endl;
//            auto run_number=find_it-heads.begin();
//            copy_space[cur_it-begin]=-1*(run_number+1); //record the run number where the current index should reside in
//            *find_it=*cur_it++;
//            run_sizes[run_number]++;
//            head_offsets[run_number]++;
////            if(find_it!=heads.begin()){
////                while(cur_it< end && *cur_it<*find_it && *cur_it>*(find_it-1)){
////                    copy_space[cur_it-begin]=-1*(run_number+1);
////                    *find_it=*cur_it++;
////                    run_sizes[run_number]++;
////                    head_offsets[run_number]++;
////
////                }
////            }
////            else{
////                while(cur_it<end && *cur_it<*find_it){
////                    copy_space[cur_it-begin]=run_number+1;
////                    *find_it=*cur_it++;
////                    run_sizes[run_number]++;
////                    head_offsets[run_number]++;
////
////                }
////            }
//            if(find_it==heads.end()-1)
//                cur_max=*(cur_it-1);
//        }
        else{
           // std::cout << "Making new run " << *cur_it << std::endl;;
            copy_space[cur_it-begin]=tails.size()+1;
            heads.push_back(*cur_it);
            tails.push_back(*cur_it++);
            head_offsets.push_back(0);
            run_sizes.push_back(1);
            //search_region=tails.size();
            search_region=std::min(tails.size(),(size_t)1000);
            cur_min=*(cur_it-1);
            cur_max=*(cur_it-1);

        }





    }
    //std::cout << "Creating runs " << boost::timer::format(creating_runs.elapsed()) << " there are " << run_sizes.size() << std::endl;
    //std::cout << "Finding runs " << boost::timer::format(finding_timer.elapsed()) << std::endl;

//    for(auto it=begin;it<end;++it)
//        std::cout << *it << " ";
//    std::cout << std::endl;

//    for(auto&x: run_sizes)
//        std::cout << x << " ";
//    std::cout << std::endl;

    //sort by run_sizes, and update the runs designations array to reflect hte new sort order
    //boost::timer::cpu_timer sorting_runs;
    SortRunsBySize(run_sizes,head_offsets,copy_space);
    //std::cout << "Sorting runs " << boost::timer::format(sorting_runs.elapsed()) << std::endl;


//    std::cout << "Sorted" << std::endl;
//    for(auto&x: run_sizes)
//        std::cout << x << " ";
//    std::cout << std::endl;
//
//    for(auto&x: copy_space)
//        std::cout << x << " ";
//    std::cout << std::endl;
    //calculate run_starts
    run_starts.resize(run_sizes.size());
    run_starts[0]=0;
    for(auto _idx=1;_idx<run_sizes.size();++_idx){
        run_starts[_idx]=run_starts[_idx-1]+run_sizes[_idx-1];

    }



    std::vector<index_type> head_starts=head_offsets;

    cur_it=begin;
    auto follow_it=followBegin;

//    for(auto it=begin;it<end;++it){
//        std::cout << *it << " " ;
//    }
//    std::cout << std::endl;
//    for(auto it=copy_space.begin();it<copy_space.end();++it){
//        std::cout << *it << " " ;
//    }
//    std::cout << std::endl;
    //pack the runs in, with smallest runs going in first
    //boost::timer::cpu_timer packing_runs;
    int cur_run_start,cur_head_offset;
    for(auto _idx=0;_idx<end-begin;++_idx,++cur_it,++follow_it){
       //std::cout << "Taking " << _idx << "run " << copy_space[_idx] << "put into " << temp_starts[copy_space[_idx]] << std::endl;
        auto cur_run=std::abs(copy_space[_idx])-1;
        if(copy_space[_idx]>0){
            cur_run_start=run_starts[cur_run]+head_offsets[cur_run];
            head_offsets[cur_run]++;
        }
        else{
            cur_run_start=run_starts[cur_run]+head_starts[cur_run]-1;
            head_starts[cur_run]--;

        }
        mainScratch[cur_run_start]=*cur_it;
        followScratch[cur_run_start]=*follow_it;
    }
    //std::cout << "Packing runs " << boost::timer::format(packing_runs.elapsed()) << std::endl;
//    for(auto it=mainScratch.begin();it<mainScratch.end();++it){
//        std::cout << *it << std::endl;
//    }
//    std::cout << "main scratch done" << std::endl;
    //boost::timer::cpu_timer merging_runs;
    int last_array=PingPongMerge(mainScratch.begin(),followScratch.begin(),begin, followBegin, run_starts, run_sizes);

    if (last_array==1){
        std::copy(mainScratch.begin(), mainScratch.end(),begin);
        std::copy(followScratch.begin(), followScratch.end(),followBegin);
    }
    //gfx::timsort(begin,end,followBegin,followBegin+(end-begin),mainScratch,followScratch);


    //std::cout << "Merging runs " << boost::timer::format(packing_runs.elapsed()) << std::endl;
    //std::cout << "Total " << boost::timer::format(total.elapsed()) << std::endl;
//    auto it2=mainScratch.begin();
//    for(auto it=begin;it<end;++it){
//        std::cout << *it << std::endl;
//    }

    //std::cout << "Finished merge " << *begin << std::endl;

}

template<typename index_type>
void SortRunsBySize(std::vector<int>& run_sizes,std::vector<index_type>& run_designations){

    std::vector<int> run_refs(run_sizes.size());
    std::vector<int> run_refs2(run_sizes.size());
    typedef decltype(run_refs.begin()) it_type;


    for(auto _idx=0;_idx<run_refs.size();++_idx){
        run_refs[_idx]=_idx;
    }
//    std::cout << "Run sizes" << std::endl;
//    for(auto &x: run_sizes)
//        std::cout << x << " ";
//    std::cout << std::endl;


    //sort runs by size
    internal::Introsort(run_sizes.begin(),run_sizes.end(),std::less<index_type>(),internal::DualSwapper<it_type,it_type>(run_sizes.begin(),run_refs.begin()));
    //reverse the run specifications
    for(auto _idx=0;_idx<run_refs.size();++_idx){
        run_refs2[run_refs[_idx]]=_idx; //run designations are 1-indexed
    }
//    std::cout << "Run refs" << std::endl;
//    for(auto &x: run_refs2)
//        std::cout << x << " ";
//    std::cout << std::endl;
    //update the run designations to the sorted layout


    for(auto it=run_designations.begin();it<run_designations.end();++it){

        *it=run_refs2[*it];

    }



}


template <typename ForwardIt,typename FollowIt>
void PingPongSortRuns(ForwardIt begin, ForwardIt end,FollowIt followBegin, std::vector<typename std::iterator_traits<ForwardIt>::value_type> & mainScratch,
                      std::vector<typename std::iterator_traits<FollowIt>::value_type> & followScratch) {


//        for(auto it=begin;it<end;++it){
//        std::cout << *it << std::endl;
//    }


    std::vector<int> run_sizes;
    run_sizes.reserve(std::floor(std::sqrt(end-begin)));
    std::vector<int> run_designations(end-begin);

    run_sizes.push_back(1);
    for(auto cur_it=begin+1;cur_it<end;++cur_it){
        if(*cur_it>=*(cur_it-1)){
            run_sizes.back()++;
            run_designations[cur_it-begin]=run_sizes.size()-1;
        }
        else{
            run_sizes.push_back(1);
            run_designations[cur_it-begin]=run_sizes.size()-1;
        }
    }
//    for(auto &x:run_sizes){
//        std::cout << x<< std::endl;
//    }
    SortRunsBySize(run_sizes,run_designations);

    std::vector<int> run_starts(run_sizes.size());
    run_starts[0]=0;
    for(auto _idx=1;_idx<run_sizes.size();++_idx){
        run_starts[_idx]=run_starts[_idx-1]+run_sizes[_idx-1];

    }
    auto temp_starts=run_starts;
    auto cur_it=begin;
    auto follow_it=followBegin;
    for(auto _idx=0;_idx<end-begin;++_idx,++cur_it,++follow_it){
       //std::cout << "Taking " << _idx << "run " << copy_space[_idx] << "put into " << temp_starts[copy_space[_idx]] << std::endl;
        auto cur_run=run_designations[_idx];
        auto cur_run_start=temp_starts[cur_run];
        mainScratch[cur_run_start]=*cur_it;
        followScratch[cur_run_start]=*follow_it;
        temp_starts[cur_run]++;
    }

    int last_array=PingPongMerge(mainScratch.begin(),followScratch.begin(),begin, followBegin, run_starts, run_sizes);
    if (last_array==1){
        std::copy(mainScratch.begin(), mainScratch.end(),begin);
        std::copy(followScratch.begin(), followScratch.end(),followBegin);
    }

//        auto it2=mainScratch.begin();
//    for(auto it=begin;it<end;++it){
//        std::cout << *it << std::endl;
//    }

}





template<typename MainIterator,typename FollowIterator,typename ScratchIterator,typename ScratchIterator2,typename HistIt,typename BufferIt, typename FollowBufferIt,unsigned long Radix,unsigned long bufferSize>
void RadixIteration(MainIterator begin, FollowIterator followBegin,ScratchIterator scratchIt1, ScratchIterator2 scratchIt2,HistIt histIt,BufferIt bufferIt,FollowBufferIt followBufferIt,size_t length,unsigned long shift_amount){

    //std::cout << "RadixIteration bitMask " << bitMask << " shift_amount " << shift_amount << std::endl;

    constexpr unsigned long bitMask=Radix-1;

    std::array<unsigned int, Radix> bufferCounts;
    for(auto &i: bufferCounts)
        i=0;


    auto followIt=followBegin;
	for (auto it=begin;it<begin+length;++it,++followIt) {



		auto pos = *it>>shift_amount & bitMask;
        //std::cout << "*it "  << *it << " pos " << pos << " histIt[pos] " << histIt[pos] << std::endl;
        *(bufferIt+pos*bufferSize+bufferCounts[pos])=*it;
        *(followBufferIt+pos*bufferSize+bufferCounts[pos]++)=*followIt;
        if(bufferCounts[pos]==bufferSize){
            std::copy(bufferIt+pos*bufferSize,bufferIt+(pos+1)*bufferSize,scratchIt1 + histIt[pos]);
            std::copy(followBufferIt+pos*bufferSize,followBufferIt+(pos+1)*bufferSize,scratchIt2 + histIt[pos]);
            histIt[pos]+=bufferSize; //update the histogram counter
            bufferCounts[pos]=0; //reset the buffer ctr
        }

	}

	for(int i=0;i<Radix;++i){
        if(bufferCounts[i]){
            std::copy(bufferIt+i*bufferSize,bufferIt+i*bufferSize+bufferCounts[i],scratchIt1 + histIt[i]);
            std::copy(followBufferIt+i*bufferSize,followBufferIt+i*bufferSize+bufferCounts[i],scratchIt2 + histIt[i]);
        }

	}

}

template<typename MainIterator,typename FollowIterator,typename ScratchIterator,typename ScratchIterator2,typename HistIt,unsigned long Radix>
void RadixIteration(MainIterator begin, FollowIterator followBegin,ScratchIterator scratchIt1, ScratchIterator2 scratchIt2,HistIt histIt,size_t length,unsigned long shift_amount){

    //std::cout << "RadixIteration bitMask " << bitMask << " shift_amount " << shift_amount << std::endl;

    constexpr unsigned long bitMask=Radix-1;


    auto followIt=followBegin;
	for (auto it=begin;it<begin+length;++it,++followIt) {



		auto pos = *it>>shift_amount & bitMask;
        //std::cout << "*it "  << *it << " pos " << pos << " histIt[pos] " << histIt[pos]+1 << std::endl;

		*(scratchIt1 + histIt[pos]) = *it;

		//std::cout << "Assigned scratchIt1" << std::endl;
		*(scratchIt2 + histIt[pos]++) = *followIt;
		//std::cout << "Assigned scratchIt2" << std::endl;
	}

}

template<typename MainIterator,typename FollowIterator,typename ScratchIterator,typename ScratchIterator2>
void RadixSort(MainIterator begin, MainIterator end, FollowIterator followBegin,ScratchIterator scratchIt1, ScratchIterator2 scratchIt2, typename std::iterator_traits<MainIterator>::value_type max_size) {
    typedef typename std::iterator_traits<MainIterator>::value_type index_type;
    typedef typename std::iterator_traits<FollowIterator>::value_type follow_type;
    static_assert(std::is_integral<index_type>::value,"Must input an integral type for main interator");




	constexpr unsigned long Radix = 2048;
	constexpr unsigned long bitMask = Radix-1;
	constexpr unsigned long bit_size = 11;
	constexpr unsigned long bufferSize=8;

	auto current_bit_size=(unsigned long)std::ceil(std::log2(max_size));
	unsigned long pass_number=std::ceil((double)current_bit_size/bit_size);

	typedef std::vector<size_t> Hist;
	Hist histograms(Radix*pass_number,0);
	typedef std::array<index_type,Radix*bufferSize> Buffer;
	Buffer buffer;
	typedef std::array<follow_type,Radix*bufferSize> FollowBuffer;
	FollowBuffer followBuffer;
	typedef typename Hist::iterator HistIt;
	typedef typename Buffer::iterator BufferIt;
	typedef typename FollowBuffer::iterator FollowBufferIt;
    //std::cout << "Histcount " << histCount << std::endl;
	auto hist_begin=histograms.begin();



    // 1.  parallel histogramming pass
	//

	for (auto it=begin; it < end; ++it) {
	    auto val=*it;
		for(int i=0;i<pass_number;++i){
            auto curHistBegin=hist_begin+i*Radix;
            size_t pos=val>>(i*bit_size) & bitMask;
            curHistBegin[pos]++;
		}

	}
	std::vector<size_t> sums(pass_number,0);



    size_t tsum;
    for(int j=0;j<pass_number;++j){
        auto curHistBegin=hist_begin+j*Radix;
        for (size_t i = 0; i < Radix; ++i) {
            tsum=*(curHistBegin+i)+sums[j];
            *(curHistBegin+i)=sums[j];
            sums[j]=tsum;
        }

    }

    size_t length=end-begin;

    int starting_pass=1; //pass number starts at one, skipping hte most LSD


    unsigned long shift_amount=starting_pass*bit_size;
    //std::cout << "pass_number " << pass_number << std::endl;
    for(int i=starting_pass;i<pass_number;++i){
//        for(auto it=histograms.begin()+i*kHist;it<histograms.begin()+(i+1)*kHist;++it)
//            std::cout << *it << std::endl;

        if ((i-starting_pass)%2==0){
            RadixIteration<MainIterator,FollowIterator,ScratchIterator,ScratchIterator2,HistIt,BufferIt,FollowBufferIt,Radix,bufferSize>(begin,followBegin,scratchIt1,scratchIt2,
                                                         histograms.begin()+i*Radix,buffer.begin(),followBuffer.begin(),length,shift_amount);
        }
        else{
            RadixIteration<MainIterator,FollowIterator,ScratchIterator,ScratchIterator2,HistIt,BufferIt,FollowBufferIt,Radix,bufferSize>(scratchIt1,scratchIt2,begin,followBegin,
                                                         histograms.begin()+i*Radix,buffer.begin(),followBuffer.begin(),length,shift_amount);
        }


       shift_amount+=bit_size;
       //std::cout << "Pass " << i << std::endl;

    }
    //if the elements remain in the scratch arrays, copy them back to the appropriate location
    if(std::max((int)pass_number-starting_pass,0)%2==1){
        std::copy(scratchIt1,scratchIt1+length,begin);
        std::copy(scratchIt2,scratchIt2+length,followBegin);
    }



    InsertionSortImproved(begin,begin+length,followBegin,std::less<index_type>()); //clean up


}

} // namespace internal

}// namespace LibMIA

#endif // LIBMIAALGORITHM
