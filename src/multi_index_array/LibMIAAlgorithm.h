// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.

/***************************************************************************
   The code for introsort and natural mergesort is modified from Keith Schwarz's versions.
   Original code can be found online at:
     http://www.keithschwarz.com/interesting/code/?dir=introsort
     http://www.keithschwarz.com/interesting/code/?dir=natural-mergesort

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

#ifndef LIBMIAALGORITHM
#define LIBMIAALGORITHM
//#include <boost/timer/timer.hpp>

namespace LibMIA
{







namespace internal
{

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
    const size_t kBlockSize = 24;


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
void InsertionSort(RandomIterator begin, RandomIterator end,
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
  Introsort(begin,end,comp,std::swap<typename std::iterator_traits<RandomIterator>::value_type>());
}

/* Non-comparator version calls the comparator version. */
template <typename RandomIterator>
void Introsort(RandomIterator begin, RandomIterator end) {
  Introsort(begin, end,
            std::less<typename std::iterator_traits<RandomIterator>::value_type>());
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



    const size_t kBlockSize = 24;
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





} // namespace internal

}// namespace LibMIA

#endif // LIBMIAALGORITHM
