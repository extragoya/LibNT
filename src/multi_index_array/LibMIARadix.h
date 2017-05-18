
#ifndef LIBMIA_RADIX_H
#define LIBMIA_RADIX_H
#include <array>
#include <vector>
#include <memory>
#include "LibMIAUtil.h"
#include "LibMIAAlgorithm.h"
#include "IndexUtil.h"
#include "libdivide.h"

namespace LibMIA
{







namespace internal
{

// Copyright(c), Victor J. Duvanenko, 2009

// Listing 3
template<unsigned int BitMask, class _Type1>
inline typename std::make_unsigned<_Type1>::type extractDigit(const  _Type1 a, const unsigned int shiftRightAmount )
{
    typedef typename std::make_unsigned<_Type1>::type unsigned_Type;
    unsigned_Type digit = (unsigned_Type)((a>>shiftRightAmount) & BitMask );
    return digit;
    // extract the digit we are sorting based on return digit;
}

template<class _Type1>
inline typename std::make_unsigned<_Type1>::type extractDigit( _Type1 a, unsigned int shiftRightAmount,unsigned int BitMask )
{
    typedef typename std::make_unsigned<_Type1>::type unsigned_Type;
    unsigned_Type digit = (unsigned_Type)((a>>shiftRightAmount) & BitMask );
    return digit;
    // extract the digit we are sorting based on return digit;
}







// Adapted from Victor J. Duvanenko's article on In-Place Radix Sort http://www.drdobbs.com/parallel/parallel-in-place-radix-sort-simplified/229000734?pgno=1

// Function template In-place MSD Radix Sort implementation (simplified).
template<unsigned long PowerOfTwoRadix ,unsigned long Log2ofPowerOfTwoRadix, long Threshold, long InsertThreshold, class ForwardIt, class FollowIt >
inline void _RadixSort_Unsigned_PowerOf2Radix_L1( ForwardIt begin, FollowIt followBegin, size_t length, unsigned long shiftRightAmount )
{

	//std::cout << "shiftRightAmount " << shiftRightAmount << std::endl;
    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;

	libmia_constexpr unsigned long BitMask = PowerOfTwoRadix - 1; //use constexpr workaround
    std::array<unsigned long,PowerOfTwoRadix> count;

    for(auto &e:count)
        e=0;

    for ( auto it=begin;it< begin+ length; ++it)
        count[((*it)>>shiftRightAmount) & BitMask ]++; // Scan the array and count the number of times each value appears

    std::array<long,PowerOfTwoRadix+1> startOfBin;
    std::array< long,PowerOfTwoRadix+1> endOfBin;
    long nextBin = 1;

    startOfBin[ 0 ] = endOfBin[ 0 ] = 0;    startOfBin[ PowerOfTwoRadix ] = -1;         // sentinal
    for( unsigned long i = 1; i < PowerOfTwoRadix; ++i )
        startOfBin[ i ] = endOfBin[ i ] = startOfBin[ i - 1 ] + count[ i - 1 ];

    for ( long _current = 0; _current < length; )
    {
        unsigned long digit;
        auto _current_element = *(begin+_current); // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
        auto _current_element_follow = *(followBegin+_current); // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
        while( endOfBin[ digit = (unsigned long)((_current_element)>>shiftRightAmount& BitMask)] != _current ){
            std::swap( _current_element, *(begin+ endOfBin[ digit ] ) );
            std::swap( _current_element_follow, *(followBegin+ endOfBin[ digit ]++ ) );
        }
        *(begin+_current) = _current_element;
        *(followBegin+_current)= _current_element_follow;

        endOfBin[ digit ]++;
        while( endOfBin[ nextBin - 1 ] == startOfBin[ nextBin ] )
            ++nextBin;   // skip over empty and full bins, when the end of the current bin reaches the start of the next bin
        _current = endOfBin[ nextBin - 1 ];
    }
    if(shiftRightAmount==0){
        return;

    }else{
        shiftRightAmount -= Log2ofPowerOfTwoRadix; //should always be a multiple of Log2ofPowerOfTwoRadix
    }



    for( unsigned long i = 0; i < PowerOfTwoRadix; ++i )
    {
        long numberOfElements = endOfBin[ i ] - startOfBin[ i ];
        if ( numberOfElements >= Threshold)     // endOfBin actually points to one beyond the bin
			_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold >(begin + startOfBin[i], followBegin + startOfBin[i], numberOfElements, shiftRightAmount);
		else if (numberOfElements >= InsertThreshold)
            introsort_detail::IntrosortRec(begin+ startOfBin[ i ],begin+ startOfBin[ i ]+numberOfElements,introsort_detail::IntrosortDepth(begin+ startOfBin[ i ],begin+ startOfBin[ i ]+numberOfElements),
                    std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin+ startOfBin[ i ],followBegin+ startOfBin[ i ]));

    }

}


// Adapted from Victor J. Duvanenko's article on In-Place Radix Sort http://www.drdobbs.com/parallel/parallel-in-place-radix-sort-simplified/229000734?pgno=1

// Function template In-place MSD Radix Sort implementation (simplified). This is one used to sort CO format, only implemented for testing and benchmarking purposes to compare against the LCO (it assumes each value in the CO format is less than 2^11)
template<unsigned long PowerOfTwoRadix, unsigned long Log2ofPowerOfTwoRadix, long Threshold, class Swapper, class ForwardIt, class FollowIt,typename LinClass,typename Comp >
inline void _RadixSort_Unsigned_PowerOf2Radix_L1(ForwardIt begin, FollowIt followBegin, size_t length,  unsigned long shiftRightAmount, const LinClass & linIdxMaker, const Comp & comp)
{

	//std::cout << "shiftRightAmount " << shiftRightAmount << std::endl;
	

	libmia_constexpr unsigned long BitMask = PowerOfTwoRadix - 1; //use constexpr workaround
	std::array<unsigned long, PowerOfTwoRadix> count;

	for (auto &e : count)
		e = 0;
	decltype (linIdxMaker(begin)) maxy=0;
	for (auto it = begin; it < begin + length; ++it){
		maxy = std::max(maxy, (linIdxMaker(it) >> shiftRightAmount) & BitMask);
		count[(linIdxMaker(it) >> shiftRightAmount) & BitMask]++; // Scan the array and count the number of times each value appears
	}

	//std::cout << "Counted " << length << " max " << maxy << " " << count[maxy] << std::endl;
	
	auto swapper = Swapper(begin, followBegin);

	std::array<long, PowerOfTwoRadix + 1> startOfBin;
	std::array< long, PowerOfTwoRadix + 1> endOfBin;
	long nextBin = 1;

	startOfBin[0] = endOfBin[0] = 0;    startOfBin[PowerOfTwoRadix] = -1;         // sentinal
	for (unsigned long i = 1; i < PowerOfTwoRadix; ++i)
		startOfBin[i] = endOfBin[i] = startOfBin[i - 1] + count[i - 1];

	//std::cout << endOfBin[maxy] << std::endl;
	auto temp_it = begin + 2;
	
	for (long _current = 0; _current < length;)
	{
		unsigned long digit;
		auto _current_element_it = (begin + _current); 
		//auto _current_element_follow = *(followBegin + _current); // get the compiler to recognize that a register can be used for the loop instead of a[_current] memory location
		//std::cout << "About to cycle length: " << length << " _current: " << _current << " digit: " << (unsigned long)((linIdxMaker(_current_element_it) >> shiftRightAmount)& BitMask) << " endOfBin[digit]: " << endOfBin[digit = (unsigned long)(linIdxMaker(_current_element_it) >> shiftRightAmount& BitMask)] << std::endl;
		//std::cout << linIdxMaker(_current_element_it) << std::endl;
		while (endOfBin[digit = (unsigned long)((linIdxMaker(_current_element_it) >> shiftRightAmount)& BitMask)] != _current){			

			/*if (endOfBin[digit] >= length){
				std::cout << "Crap!" << std::endl;
				std::cout << " digit: " << digit << " endOfBin[digit]: " << endOfBin[digit] << " linIdxMaker: " << linIdxMaker(_current_element_it) << std::endl;
			}*/
			auto it = begin + endOfBin[digit];
			auto old = linIdxMaker(_current_element_it);
			auto old2 = linIdxMaker(it);
			//std::cout << *(_current_element_it.m_iter) << " " << *(_current_element_it.m_iter + 1) << " " << *(it.m_iter) << " " << *(it.m_iter + 1) << std::endl;
			
			swapper(_current_element_it, (begin + endOfBin[digit])); //will also swap follower
			//std::cout << *(_current_element_it.m_iter) << " " << *(_current_element_it.m_iter + 1) << " " << *(it.m_iter) << " " << *(it.m_iter + 1) << std::endl;
			

			endOfBin[digit]++;
			//std::cout << *(_current_element_it.m_iter) << " " << *(_current_element_it.m_iter + 1) << " " << *(it.m_iter) << " " << *(it.m_iter + 1) << std::endl;
			//exit(0);
			//std::swap(_current_element_follow, *(followBegin + endOfBin[digit]++));
		}
		
		//*(followBegin + _current) = _current_element_follow;
		
		endOfBin[digit]++;
		while (endOfBin[nextBin - 1] == startOfBin[nextBin])
			++nextBin;   // skip over empty and full bins, when the end of the current bin reaches the start of the next bin
		_current = endOfBin[nextBin - 1];
		//std::cout << "Finished cylcing " << nextBin << std::endl;
	}
	//std::cout << "Performed radix pass" << std::endl;
	if (shiftRightAmount == 0){
		return;

	}
	else{
		shiftRightAmount -= Log2ofPowerOfTwoRadix; //should always be a multiple of Log2ofPowerOfTwoRadix
	}



	for (unsigned long i = 0; i < PowerOfTwoRadix; ++i)
	{
		long numberOfElements = endOfBin[i] - startOfBin[i];
		if (numberOfElements >= Threshold)     // endOfBin actually points to one beyond the bin
			_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold,  Swapper > (begin + startOfBin[i], followBegin + startOfBin[i], numberOfElements, shiftRightAmount, linIdxMaker, comp);
		else if (numberOfElements >= 50)
			introsort_detail::IntrosortRec(begin + startOfBin[i], begin + startOfBin[i] + numberOfElements, introsort_detail::IntrosortDepth(begin + startOfBin[i], begin + startOfBin[i] + numberOfElements),
			comp, Swapper(begin + startOfBin[i], followBegin + startOfBin[i]));

	}

}

//MSD Radix Sort, this version is used to text the sorting for CO format. Implemented to test and graph sorting performance between CO and LCO formats
template< class ForwardIt, class FollowIt,typename linIndexType,size_t co_length>
inline void RadixSortInPlace_PowerOf2Radix_Unsigned(ForwardIt  begin, FollowIt  followBegin, size_t length, size_t max_size, const std::array<size_t, co_length> & sortOrder, const std::array<linIndexType, co_length>& dims)
{

	/*using namespace std::chrono;
	typedef std::chrono::duration<float> float_seconds;
	high_resolution_clock::time_point t1, t2;*/
	//t1 = high_resolution_clock::now();
	typedef typename std::iterator_traits<ForwardIt>::value_type index_type;

    const long Threshold                      =  3000;
    const unsigned long PowerOfTwoRadix       = 2048;    // Radix - must be a power of 2
    const unsigned long Log2ofPowerOfTwoRadix =   11;    // log( 2048 ) = 11
    auto current_bit_size=(unsigned long)std::ceil(std::log2(max_size));

	auto number_shifts = (unsigned long)std::ceil((double)current_bit_size / Log2ofPowerOfTwoRadix);

	unsigned long shiftRightAmount;
	if (current_bit_size >= Log2ofPowerOfTwoRadix)
		shiftRightAmount = (number_shifts - 1)*Log2ofPowerOfTwoRadix;
	else
		shiftRightAmount = 0;

  
	auto linIdxMaker = [sortOrder, dims](ForwardIt master_it){
		auto it = master_it.m_iter;
		linIndexType ret = *(it+sortOrder[0]);
		auto mult = dims[sortOrder[0]];
		static_for<1,co_length>::_do([&](int i)
		{
			ret += *(it + sortOrder[i])*mult;
			mult *= dims[sortOrder[i]];
		});
		return ret;
	};

	


	
	auto sort_compare = [sortOrder](const index_type& idx1, const index_type& idx2){
		const index_type* it = &idx1;
		const index_type* it2 = &idx2;
		for (int i = co_length-1; i >= 0; --i){
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

	
	typedef CoSwapper<ForwardIt, FollowIt, co_length> Swapper;
	
    if ( length >= Threshold )
		_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold,Swapper >(begin, followBegin, length, shiftRightAmount , linIdxMaker, sort_compare);
    else
        introsort_detail::IntrosortRec(begin,begin+length,introsort_detail::IntrosortDepth(begin,begin+length),
		sort_compare, Swapper(begin,followBegin));
	
	
	//t2 = high_resolution_clock::now();
	//std::cout << "Straight Before Insert:"  << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
	//t1 = high_resolution_clock::now();
    //clean up the last
    //InsertionSort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin,followBegin));
	InsertionSort(begin, begin + length,  sort_compare, Swapper(begin, followBegin));
	//t2 = high_resolution_clock::now();
	//std::cout << "Straight Insert:" << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
}


//MSD Radix Sort, this version is used to text the sorting for CO format, for the MTT version of the CO format. Implemented to test and graph sorting performance between CO and LCO formats
template< class ForwardIt, class FollowIt, typename linIndexType, size_t co_length>
inline void RadixSortInPlace_PowerOf2Radix_Unsigned_MTT(ForwardIt  begin, FollowIt  followBegin, size_t length, size_t max_size, const std::array<size_t, co_length> & sortOrder, const std::array<linIndexType, co_length>& dims)
{

	/*using namespace std::chrono;
	typedef std::chrono::duration<float> float_seconds;
	high_resolution_clock::time_point t1, t2;*/
	//t1 = high_resolution_clock::now();
	

	const long Threshold = 3000;
	const unsigned long PowerOfTwoRadix = 2048;    // Radix - must be a power of 2
	const unsigned long Log2ofPowerOfTwoRadix = 11;    // log( 2048 ) = 11
	auto current_bit_size = (unsigned long)std::ceil(std::log2(max_size));

	auto number_shifts = (unsigned long)std::ceil((double)current_bit_size / Log2ofPowerOfTwoRadix);

	unsigned long shiftRightAmount;
	if (current_bit_size >= Log2ofPowerOfTwoRadix)
		shiftRightAmount = (number_shifts - 1)*Log2ofPowerOfTwoRadix;
	else
		shiftRightAmount = 0;


	


	auto linIdxMaker = [sortOrder, dims](ForwardIt master_it){
		auto it_array = master_it.m_iter_array;
		linIndexType ret = *(it_array[sortOrder[0]]);
		auto mult = dims[sortOrder[0]];
		static_for<1, co_length>::_do([&](int i)
		{

			ret += *(it_array[sortOrder[i]])*mult;
			mult *= dims[sortOrder[i]];


		});
		
		return ret;
	};


	typedef typename std::iterator_traits<ForwardIt>::value_type it_value_type;
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


	typedef CoSwapperMTT<ForwardIt, FollowIt, co_length> Swapper;

	if (length >= Threshold)
		_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, Swapper >(begin, followBegin, length, shiftRightAmount, linIdxMaker, sort_compare);
	else
		introsort_detail::IntrosortRec(begin, begin + length, introsort_detail::IntrosortDepth(begin, begin + length),
		sort_compare, Swapper(begin, followBegin));


	//t2 = high_resolution_clock::now();
	//std::cout << "Straight Before Insert:"  << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
	//t1 = high_resolution_clock::now();
	//clean up the last
	//InsertionSort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin,followBegin));
	InsertionSort(begin, begin + length, sort_compare, Swapper(begin, followBegin));
	//t2 = high_resolution_clock::now();
	//std::cout << "Straight Insert:" << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
}

//MSD Radix Sort
template< class ForwardIt, class FollowIt>
inline void RadixSortInPlace_PowerOf2Radix_Unsigned(ForwardIt  begin, FollowIt  followBegin, size_t length, size_t max_size)
{

	/*using namespace std::chrono;
	typedef std::chrono::duration<float> float_seconds;
	high_resolution_clock::time_point t1, t2;*/
	//t1 = high_resolution_clock::now();
	typedef typename std::iterator_traits<ForwardIt>::value_type index_type;

	const long Threshold = 3000;
	const long InsertThreshold = 50;
	const unsigned long PowerOfTwoRadix = 2048;    // Radix - must be a power of 2
	const unsigned long Log2ofPowerOfTwoRadix = 11;    // log( 2048 ) = 11
	auto current_bit_size = (unsigned long)std::ceil(std::log2(max_size));

	auto number_shifts = (unsigned long)std::ceil((double)current_bit_size / Log2ofPowerOfTwoRadix);

	unsigned long shiftRightAmount;
	if (current_bit_size >= Log2ofPowerOfTwoRadix)
		shiftRightAmount = (number_shifts - 1)*Log2ofPowerOfTwoRadix;
	else
		shiftRightAmount = 0;


	if (length >= Threshold)
		_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold >(begin, followBegin, length, shiftRightAmount);
	else
		introsort_detail::IntrosortRec(begin, begin + length, introsort_detail::IntrosortDepth(begin, begin + length),
		std::less<index_type>(), internal::DualSwapper<ForwardIt, FollowIt>(begin, followBegin));


	//t2 = high_resolution_clock::now();
	//std::cout << "Straight Before Insert:"  << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
	//t1 = high_resolution_clock::now();
	//clean up the last
	//InsertionSort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin,followBegin));
	InsertionSortImproved(begin, begin + length, followBegin, std::less<index_type>());
	//t2 = high_resolution_clock::now();
	//std::cout << "Straight Insert:" << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
}



template< unsigned int PowerOfTwoRadix,class ForwardIt,class FollowIt,class ScratchIt,class ScratchItFollow>
inline std::array<size_t,PowerOfTwoRadix+1> _performRadixExtractDigit( ForwardIt  begin1, FollowIt  followBegin1,ScratchIt  begin2, ScratchItFollow  followBegin2, size_t length,
                        unsigned long cur_shiftRightAmount ){



    constexpr unsigned int numberOfBins=PowerOfTwoRadix;

    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;

    typedef typename std::make_unsigned<index_type>::type unsigned_Type;

    constexpr unsigned int BitMask =PowerOfTwoRadix-1;


    std::array<size_t,PowerOfTwoRadix+1> count;


    for(auto &i:count)
        i=0;


    for ( auto it=begin1;it<begin1+length;++it ) // Scan the array and count the number of times each value in the current bitMask appears
        count[ extractDigit<BitMask>( ((unsigned_Type)*it), cur_shiftRightAmount )+1 ]++;

    for(auto idx=1;idx<numberOfBins+1;++idx) // Turn the counts into starting offsets into the array
        count[idx]+=count[idx-1];

    auto followIt=followBegin1;
    auto temp_count=count;

    for ( auto it=begin1;it<begin1+length;++it,++followIt){ //perform radix sort on the current bitMask, using the starting offsets
        auto cur_Digit=extractDigit<BitMask>( ((unsigned_Type)*it), cur_shiftRightAmount );
        *(begin2+temp_count[cur_Digit])=*it;
        *(followBegin2+temp_count[cur_Digit]++)=*followIt;


    }
    return count;


}








//!Adpated from  Victor J. Duvanenko's Stable MSB Radix source code found at http://www.drdobbs.com/tools/algorithm-improvement-through-performanc/222200161
// Altered to suit the needs of sparse MIA permutations
//
template< unsigned long PowerOfTwoRadix, unsigned long Log2ofPowerOfTwoRadix, long Threshold,class ForwardIt,class FollowIt,class ScratchIt,class ScratchItFollow>
inline void _RadixSort_Sort_Simple( ForwardIt  begin1, FollowIt  followBegin1,ScratchIt  begin2, ScratchItFollow  followBegin2, size_t length,
                        unsigned long cur_shiftRightAmount,bool inputArrayIsDestination){


    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;



    //std::vector<size_t> count;

    bool stageDone=false;


    auto count=_performRadixExtractDigit<PowerOfTwoRadix>(begin1,followBegin1,begin2,followBegin2,length,cur_shiftRightAmount);


    if(cur_shiftRightAmount==0){



        if ( !inputArrayIsDestination ){

            std::copy(begin2,begin2+length,begin1);
            std::copy(followBegin2,followBegin2+length,followBegin1);
        }
        return;

    }
    else{


        cur_shiftRightAmount -= Log2ofPowerOfTwoRadix; //should always be a multiple of Log2ofPowerOfTwoRadix

    }






    inputArrayIsDestination = !inputArrayIsDestination;
    for( size_t i = 0; i < count.size()-1; ++i )
    {
        size_t offset=count[i];
        size_t numOfElements = count[ i+1 ] - count[ i ];

        if ( numOfElements >= Threshold ){
            //std::cout << "Doing a radix recursion i " << i << " " << offset << " " << numOfElements << std::endl;

            _RadixSort_Sort_Simple<PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold>( begin2+offset,followBegin2+offset,begin1+offset,followBegin1+offset,numOfElements , cur_shiftRightAmount, inputArrayIsDestination);

        }
        else {
            if (numOfElements>50) {



                introsort_detail::IntrosortRec(begin2+offset,begin2+offset+numOfElements,introsort_detail::IntrosortDepth(begin2+offset,begin2+offset+numOfElements),
                    std::less<index_type>(),internal::DualSwapper<ScratchIt,ScratchItFollow>(begin2+offset,followBegin2+offset));
            }
            if ( inputArrayIsDestination ){ //normally this would be !inputArrayIsDestination, but we inverted the value above
                std::copy(begin2+offset,begin2+offset+numOfElements,begin1+offset);
                std::copy(followBegin2+offset,followBegin2+offset+numOfElements,followBegin1+offset);
            }
        }
    }


}
template<unsigned int  Threshold,unsigned int IntroThreshold, unsigned int Radix>
inline int estimateNumberOfRuns(size_t & length, size_t max_size){

	int numberOfRuns = 0;
	int radixRuns = static_cast<int>(std::ceil(std::log2(max_size) / std::log2(Radix)));
	size_t first_run = max_size / static_cast<size_t>(std::pow(Radix, radixRuns - 1));
	//std::cout << "length " << length << "First Run " << first_run << " radixRuns " << radixRuns << " max_size " << max_size <<std::endl;
	if (length > Threshold){
		length /= first_run;
		radixRuns--;
		numberOfRuns++;
		//std::cout << "did first run length " << length << "numberOfRuns " << numberOfRuns << std::endl;
		while (length > Threshold && radixRuns > 0){
			numberOfRuns++;
			radixRuns--;
			length /= Radix;
			//std::cout << "did normal run length " << length << "numberOfRuns " << numberOfRuns << std::endl;
		}
	}
	//std::cout << "radixRuns" << radixRuns << std::endl;
	if (radixRuns > 0 && length <= Threshold){
		if (length > IntroThreshold)
			numberOfRuns += 2;
		else if (length > 0)
			numberOfRuns++;
		length = 0;
	}


	//std::cout << "Number of Runs " << numberOfRuns << std::endl;
    return numberOfRuns;
}




template<unsigned int IntroThreshold, unsigned int Radix>
inline size_t getLengthThreshold(size_t max_size){

    assert(max_size);
    size_t thresholdLength=50;
    if (max_size==1)
        return thresholdLength;
	int radixRuns = static_cast<int>(std::ceil(std::log2(max_size) / std::log2(Radix)));
    size_t first_run=max_size/static_cast<size_t>(std::pow(Radix,radixRuns-1));

    if(radixRuns>1){
        for(int i=0;i<radixRuns-2;++i){
            thresholdLength*=Radix;
        }

    }
    thresholdLength*=first_run;
    return thresholdLength;
}








template<unsigned int PowerOfTwoRadix,unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,class ForwardIt,class FollowIt,class ScratchIt,class ScratchItFollow,class unsigned_Type>
inline void RadixSortStraight( ForwardIt  begin, FollowIt  followBegin, ScratchIt  scratchBegin, ScratchItFollow  scratchBegin2,size_t length,const unsigned_Type max_size)
{

    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;


    if(length<Threshold)
        Introsort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin,followBegin));


    constexpr unsigned int bufferSize=16;

    auto current_bit_size=(unsigned int)std::ceil(std::log2(max_size));
    auto number_shifts=(unsigned int)std::ceil((double)current_bit_size/Log2ofPowerOfTwoRadix);
    unsigned int shiftRightAmount;
    if(current_bit_size>=Log2ofPowerOfTwoRadix)
        shiftRightAmount= (number_shifts-1)*Log2ofPowerOfTwoRadix;
    else
        shiftRightAmount=0;


   _RadixSort_Sort_Simple<PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold>(begin,followBegin,scratchBegin,scratchBegin2,length,shiftRightAmount,false);


    //clean up the last
    InsertionSortImproved(begin,begin+length,followBegin,std::less<index_type>());


}

template<class index_type, class data_type,unsigned int PowerOfTwoRadix=2048, unsigned int Log2ofPowerOfTwoRadix=11,unsigned int Threshold=3000,unsigned int InsertThreshold=50>
class RadixShuffle
{

private:

    typedef typename std::make_unsigned<index_type>::type unsignedType;

    const std::vector<unsignedType> mMaxDims;
    const std::vector<unsignedType> mDivisors;
    std::vector<libdivide::divider<unsignedType>> mFastDivisors;
    std::vector<unsigned int> mShiftRightAmounts;
	std::vector<unsigned int> mDivisorsShiftAmount;
	std::unique_ptr<index_type[]> mScratch1;
	std::unique_ptr<data_type[]> mScratch2;
    
    std::array<size_t,PowerOfTwoRadix+1> mCountBuffer;
    std::array<size_t,PowerOfTwoRadix+1> mCountBuffer2;
    const size_t mTotalMax;
    unsigned int mTotalBitSize;
    unsigned int mTotalRightAmount;
    const bool mFirstSortOrFind;
	size_t curBufferLength;
    typedef index_type* ScratchIt1;
    typedef data_type* ScratchIt2;
public:

    /**
   * Constructor: RadixShuffle: performs a radix sort customized to work for sparse MIA shuffles. Essentially, when performing an sparse shuffle, or a permulation of a sparseMIA's
                    linear indices (what index order it uses to calculate a linear index), the resulting data has some structure that can be exploited for sorting purposes.


        \param[in] _MaxDims vector that provides the max sizes of each stage in order to create appropriate bitmasks.
        \param[in] _Divisors vector that provides the divisors needed to extract the digits for each stage
        \param[in] _TotalMax the max size of the indices
        \param[in] _SortOrFind boolean value that indicates what the first stage is - sort or find, true or false respectively

   */
    RadixShuffle(const std::vector<unsignedType> & _MaxDims,const std::vector<unsignedType> & _Divisors,const size_t _TotalMax, const bool _SortOrFind)
		:mMaxDims(_MaxDims), mDivisors(_Divisors), mTotalMax(_TotalMax), mFirstSortOrFind(_SortOrFind), curBufferLength(0)
    {
        assert(mMaxDims.size()==mDivisors.size());
        auto number_of_stages=mMaxDims.size();
        mShiftRightAmounts.resize(mMaxDims.size());
        mFastDivisors.resize(mMaxDims.size());
		mDivisorsShiftAmount.resize(mMaxDims.size());
        mTotalBitSize=(unsigned int)std::ceil(std::log2(mTotalMax));
        auto total_number_shifts=(unsigned int)std::ceil((double)mTotalBitSize/Log2ofPowerOfTwoRadix);
        if(mTotalBitSize>=Log2ofPowerOfTwoRadix)
            mTotalRightAmount= (total_number_shifts-1)*Log2ofPowerOfTwoRadix;
        else
            mTotalRightAmount=0;
        //std::array<size_t,PowerOfTwoRadix+1> count;
        for(auto idx=0; idx<number_of_stages; ++idx)
        {
            auto current_bit_size=(unsigned int)std::ceil(std::log2(mMaxDims[idx]));
            auto number_shifts=(unsigned int)std::ceil((double)current_bit_size/Log2ofPowerOfTwoRadix);


            if(current_bit_size>=Log2ofPowerOfTwoRadix)
                mShiftRightAmounts[idx]= (number_shifts-1)*Log2ofPowerOfTwoRadix;
            else
                mShiftRightAmounts[idx]=0;

            //also create a vector of libdivide::divisors
            mFastDivisors[idx]=libdivide::divider<unsignedType>(mDivisors[idx]);
			
			mDivisorsShiftAmount[idx] = std::floor(std::log2(mDivisors[idx]));
			//std::cout << "divisor " << mDivisors[idx] << " shiftDivisor " << mDivisorsShiftAmount[idx] << std::endl;
        }
		
    }

    void mResize(size_t length){
		if (curBufferLength<length){
			mScratch1.reset(new index_type[length]);
			mScratch2.reset( new data_type[length]);
        }
    }
	
	
	inline int estimateNumberOfRunsShuffle(size_t t_length, int start_index=0){


		int num_runs = 0;
		for (int i = start_index; i < mMaxDims.size(); i += 2){
			size_t cur_max_dim = (size_t)(mMaxDims[i] * mDivisors[i]) >> mDivisorsShiftAmount[i]; //account for the fact that we use truncated divisors for the power of two trick
			//std::cout << cur_max_dim << " " << mMaxDims[i] << " " << mDivisors[i] << " " << mDivisorsShiftAmount[i] << std::endl;
			num_runs += estimateNumberOfRuns<Threshold, InsertThreshold, PowerOfTwoRadix>(t_length, cur_max_dim);
			if (!t_length)
				break;
			if (i + 1 < mMaxDims.size()){
				t_length = std::max(std::ceil((float)t_length / mMaxDims[i + 1]), float(2));
			}
		}


		//std::cout << "Number of Runs " << numberOfRuns << std::endl;
		return num_runs;
	}


    //!Perform the permutation
    template<class RandomIt, class FollowIt>
    void permute(RandomIt begin, FollowIt  followBegin,size_t length)
    {
        
		
		using namespace std::chrono;
		typedef std::chrono::duration<float> float_seconds;
		high_resolution_clock::time_point t1, t2,t3,t4;
		//t3 = high_resolution_clock::now();
		if(length<Threshold){
			//std::cout << length << " < " << Threshold << std::endl;
			Introsort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<RandomIt,FollowIt>(begin,followBegin));
            return;
        }


        if(mFirstSortOrFind==false)  //if we are currently finding and not sorting
        {

			
			
			
			/*_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold>(begin, followBegin, length, mTotalRightAmount);
			
			InsertionSortImproved(begin, begin + length, followBegin, std::less<index_type>());*/
			
			if (mMaxDims.size()>1)
				RadixShuffleFind(begin,followBegin,length);
			

        }
        else
        {
            
			
			int num_runs = estimateNumberOfRunsShuffle(length);
			size_t t_length = length;
			int num_straight_runs = estimateNumberOfRuns<Threshold, InsertThreshold, PowerOfTwoRadix>(t_length, mTotalMax);
			//std::cout << "max size " << mMaxDims.front() << " length " << length << " estimate of runs " << num_runs << std::endl;
			//std::cout << "total_max_size " << mTotalMax << " estimate of runs " << num_straight_runs << std::endl;

			if (num_runs < num_straight_runs)
            {
				//std::cout << "Doing shuffle " << std::endl;
				t1 = high_resolution_clock::now();
                mResize(length);
				t2 = high_resolution_clock::now();
				
				//std::cout << "Resize Time:" << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
				//t1 = high_resolution_clock::now();
                RadixShuffleSort(begin,followBegin,mScratch1.get(),mScratch2.get(),length,mShiftRightAmounts[0],0,false);	
				//InsertionSortImproved(begin, begin + length, followBegin, std::less<index_type>());

            }
            else
            {
                //std::cout << "Doing straight radix - not worth it " << mTotalRightAmount << std::endl;
                //RadixSortInPlace_PowerOf2Radix_Unsigned( begin,  followBegin,length,total_max_size );
				t1 = high_resolution_clock::now();
				_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>( begin,followBegin, length, mTotalRightAmount );
				t2 = high_resolution_clock::now();

				//std::cout << "Radix time Time:" << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
                //RadixSortStraight<PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold>(  begin,  followBegin, scratch1.begin(),scratch2.begin(),length,total_max_size);
				t1 = high_resolution_clock::now();
				InsertionSortImproved(begin, begin + length, followBegin, std::less<index_type>());
				t2 = high_resolution_clock::now();
				//std::cout << "Insert Time:" << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
            }

			

        }
		
		
    }
private:
    //only should run as the first recursion
    template<class RandomIt, class FollowIt>
    void RadixShuffleFind(RandomIt begin, FollowIt followBegin,size_t length);

    //only should run as later recurions
    template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
	void RadixShuffleFind(RandomIt begin, FollowIt followBegin, RandomIt2 begin2, FollowIt2 followBegin2, size_t length, int stage_index,bool inputArrayIsDestination);

	

    template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
	void RadixShuffleSort(RandomIt begin, FollowIt followBegin, RandomIt2 begin2, FollowIt2 followBegin2, size_t length, unsigned int curShiftRightAmount, int stage_index, bool inputArrayIsDestination);

    template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
    void performRadixExtractDigit(RandomIt begin1,FollowIt followBegin1,RandomIt2 begin2,FollowIt2 followBegin2,size_t length,unsigned int cur_shiftRightAmount,int stage_index);

	template<class RandomIt, class FollowIt, class RandomIt2, class FollowIt2>
	void performRadixExtractDigitRotate(RandomIt begin1, FollowIt followBegin1, RandomIt2 begin2, FollowIt2 followBegin2, const size_t length, const unsigned int cur_shiftRightAmount, const int stage_index);


    void getNextShifts(int stage_index,unsigned int & straightshiftRightAmount,unsignedType & straightEliminatorMult,unsignedType & nextEliminatorMult,size_t & shuffleLengthThreshold );

	
	

};

template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
    getNextShifts(int stage_index,unsigned int & straightshiftRightAmount,unsignedType & straightEliminatorMult,unsignedType & nextEliminatorMult,size_t & shuffleLengthThreshold )
{

    auto cur_num_divisor=mDivisors[stage_index];
    auto nextMaxSize=mMaxDims[stage_index+1]; //max size if we do divisor trick to reduce radix wordlength for next stage
    auto newStraightMaxSize=nextMaxSize*mDivisors[stage_index+1]; //max size if we just do straight radix sort for next stage

    bool isPowerOfTwo;         // if the next max sizes are a power of two, then we skip using an eliminator altogether
    straightEliminatorMult=0;
    nextEliminatorMult=0;
    isPowerOfTwo = nextMaxSize && !(nextMaxSize & (nextMaxSize - 1));
    if(!isPowerOfTwo)
        nextEliminatorMult=cur_num_divisor;
    isPowerOfTwo = newStraightMaxSize && !(newStraightMaxSize & (newStraightMaxSize - 1));
    if(!isPowerOfTwo)
        straightEliminatorMult=cur_num_divisor;
    //calculate the next shift amount if we are doing a straight radix sort

    auto next_bit_size=(unsigned int)std::ceil(std::log2(newStraightMaxSize));
    auto number_shifts=(unsigned int)std::ceil((double)next_bit_size/Log2ofPowerOfTwoRadix);

    if(next_bit_size>=Log2ofPowerOfTwoRadix)
        straightshiftRightAmount= (number_shifts-1)*Log2ofPowerOfTwoRadix;
    else
        straightshiftRightAmount=0;


    shuffleLengthThreshold= getLengthThreshold<InsertThreshold,PowerOfTwoRadix>(nextMaxSize);
}


template<class RandomIt, class FollowIt, class index_type, class data_type>
inline void SiftDown(RandomIt it, FollowIt followIt, RandomIt begin, index_type curElement, data_type followElement) {

	*it = *(it - 1);
	*followIt = *(followIt - 1);
	--it;
	--followIt;
	for (; it>begin && curElement < *(it - 1); --it, --followIt)
	{
		*it = *(it - 1);
		*followIt = *(followIt - 1);

	}
	*it = curElement;    // always necessary work/write
	*followIt = followElement;    // always necessary work/write




	

}


template<class RandomIt, class FollowIt, class RandomIt2, class FollowIt2>
inline void InsertionSortWithMove(RandomIt begin, FollowIt followIt, RandomIt2 begin2, FollowIt2 followIt2, size_t length, bool destinationIsFirst) {

	if (!length)
		return;
	
	typedef typename std::iterator_traits<RandomIt>::value_type index_type;
	if (destinationIsFirst){
		InsertionSortImproved(begin, begin + length, followIt, std::less<index_type>());
		return;
	}
	
	
	auto end = begin + length;
	auto it2 = begin2;
	*begin2 = *begin;
	*followIt2 = *followIt;
	begin++;
	it2++;
	followIt++;
	followIt2++;
	for (; begin < end; ++begin, ++followIt, ++it2, ++followIt2){
		if (*begin < *(it2 - 1))       // no need to do (j > 0) compare for the first iteration
		{
			SiftDown(it2, followIt2, begin2, *begin, *followIt);

		}
		else{
			*it2 = *begin;
			*followIt2 = *followIt;
		}


	}

}



template <typename RandomIterator, typename FollowIt, typename index_type>
inline void  InsertionSortImprovedEliminator(RandomIterator begin, RandomIterator end, FollowIt followBegin, index_type eliminator) {

	if (end==begin)
		return;
	*begin += eliminator;
	auto followIt = followBegin + 1;
	for (auto it = begin + 1; it<end; ++it, ++followIt)
	{
		*it += eliminator;
		if (*it< *(it - 1))       // no need to do (j > 0) compare for the first iteration
		{
			auto curElement = *it;
			auto followElement = *followIt;
			*it = *(it - 1);
			*followIt = *(followIt - 1);
			auto it2 = it - 1;
			auto followIt2 = followIt - 1;
			for (; it2>begin && curElement < *(it2 - 1); --it2, --followIt2)
			{
				*it2 = *(it2 - 1);
				*followIt2 = *(followIt2 - 1);

			}
			*it2 = curElement;    // always necessary work/write
			*followIt2 = followElement;    // always necessary work/write
		}
		// Perform no work at all if the first comparison fails - i.e. never assign an element to itself!
	}
}



template<class RandomIt, class FollowIt, class RandomIt2, class FollowIt2, class index_type>
inline void InsertionSortWithMoveAndEliminator(RandomIt begin, FollowIt followIt, RandomIt2 begin2, FollowIt2 followIt2, size_t length, index_type eliminator,bool destinationIsFirst) {

	if (!length)
		return;
	if (destinationIsFirst){
		InsertionSortImprovedEliminator(begin, begin + length, followIt, eliminator);
		return;
	}

	*begin += eliminator;
	auto end = begin + length;
	auto it2 = begin2;
	*begin2 = *begin;
	*followIt2 = *followIt;
	begin++;
	it2++;
	followIt++;
	followIt2++;
	for (; begin < end; ++begin, ++followIt, ++it2, ++followIt2){
		*begin += eliminator;
		if (*begin < *(it2 - 1))       // no need to do (j > 0) compare for the first iteration
		{
			SiftDown(it2, followIt2, begin2, *begin, *followIt);

		}
		else{
			*it2 = *begin;
			*followIt2 = *followIt;
		}


	}

}



template<class RandomIt,class index_type>
inline void AddEliminator(RandomIt begin, RandomIt end,  index_type eliminator){

	for (; begin < end; ++begin)
		*begin += eliminator;

}




template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
RadixShuffleFind(RandomIt begin, FollowIt followBegin, RandomIt2 begin2, FollowIt2 followBegin2, size_t length, int stage_index, bool inputArrayIsDestination)
{


    
	

	if(stage_index>=mFastDivisors.size()-1) //if this is the last stage (which is also shouldn't happen), just return, as there's no reason to do a find
        return;

    auto curIt=begin+1;
    auto cur_divisor=mFastDivisors[stage_index]; //current divisors
    auto cur_num_divisor=mDivisors[stage_index];
	index_type curValue = ((unsignedType)*begin) / cur_divisor; //get current value that defines to block to find
	
	index_type nextValue=(curValue+1)*cur_num_divisor; //get the next value that defines the end of current block

    size_t curOffset=0; //offset and number of elements in current block
    size_t curNumOfElements=1;
    index_type eliminator;
    unsigned int straightshiftRightAmount;
    unsignedType straightEliminatorMult, nextEliminatorMult;
    size_t shuffleLengthThreshold;
    getNextShifts(stage_index,straightshiftRightAmount,straightEliminatorMult,nextEliminatorMult,shuffleLengthThreshold );
    //std::cout << " nextMaxSize " << nextMaxSize << " newStraightMaxSize " << newStraightMaxSize << std::endl;
    //std::cout << "next_bit_size " << next_bit_size << " number_shifts " << number_shifts << " straightshiftRightAmount " << straightshiftRightAmount << " shuffleLengthThreshold " << shuffleLengthThreshold << std::endl;
	//std::cout << "****Cur Value***** " << curValue;
	
	
	
    for(;curIt<begin+length;++curIt){
		
        if(*curIt>=nextValue){




            if ( curNumOfElements >= Threshold ){ //if we've have enough elements to perform a radix sort

				
                if(curNumOfElements>shuffleLengthThreshold) //shuffeLengthThreshold determines how big size must be to perform Radix Permute, vs just straight radix sort
                    eliminator=curValue*nextEliminatorMult;
                else
                    eliminator=curValue*straightEliminatorMult;

                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it) //all sorting can be done modulo current index
                    *it-=eliminator;

                if(curNumOfElements>shuffleLengthThreshold){
					//perform radix permute
					//std::cout << " Suffle " << curNumOfElements << std::endl;
					RadixShuffleSort(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, mShiftRightAmounts[stage_index + 1], stage_index + 1, inputArrayIsDestination);
					AddEliminator(begin + curOffset, begin + curOffset + curNumOfElements, eliminator);

                }
                else{
					//perform straight sort
					//std::cout << " STraight " << curNumOfElements << std::endl;
                    _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>( begin+curOffset,followBegin+curOffset, curNumOfElements, straightshiftRightAmount);                    
					InsertionSortWithMoveAndEliminator(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, eliminator, !inputArrayIsDestination);					
					//if (inputArrayIsDestination)
					//InsertionSortWithMove(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, !inputArrayIsDestination);
					
                }
            }
            else{
				//std::cout << " Intro/Insert " << curNumOfElements << std::endl;
				if (curNumOfElements > InsertThreshold){ //otherwise, if we have enough to perform IntroSort, run it, otherwise do nothing, and let insertion sort do a final clean up at the end
					//std::cout << "Intro" << std::endl;
					//_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold>(begin + curOffset, followBegin + curOffset, curNumOfElements, straightshiftRightAmount);
					introsort_detail::IntrosortRec(begin + curOffset, begin + curOffset + curNumOfElements, introsort_detail::IntrosortDepth(begin + curOffset, begin + curOffset + curNumOfElements),
						std::less<index_type>(), internal::DualSwapper<RandomIt, FollowIt>(begin + curOffset, followBegin + curOffset));
				}
				/*if (inputArrayIsDestination){
					std::copy(begin + curOffset, begin + curOffset+curNumOfElements, begin2 + curOffset);
					std::copy(followBegin + curOffset, followBegin +curOffset + curNumOfElements, followBegin2 + curOffset);
				}*/
					InsertionSortWithMove(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, !inputArrayIsDestination);
					//InsertionSortImproved(begin + curOffset, begin + curOffset + curNumOfElements, followBegin + curOffset, std::less<index_type>());
            }			

            curOffset=curIt-begin;
            curNumOfElements=1;
			if (*curIt < nextValue+cur_num_divisor){ //we can avoid an integer division if we know next index value is consecutive to current one
				curValue += 1;
				nextValue += cur_num_divisor;
				
			}
			else{
				curValue = ((unsignedType)*curIt) / cur_num_divisor;// cur_divisor;
				nextValue = (curValue + 1)*cur_num_divisor;
			}
			//std::cout << "****Cur Value***** " << curValue;
			//std::cout << "cur it" << *curIt << std::endl;
        }
        else{
            ++curNumOfElements;
			//std::cout << "cur it" << *curIt << std::endl;
        }


    }
	
    //perform last section
    if ( curNumOfElements >= Threshold )
    {
        if(curNumOfElements>shuffleLengthThreshold)
            eliminator=curValue*nextEliminatorMult;
        else
            eliminator=curValue*straightEliminatorMult;

        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
            *it-=eliminator;

        if(curNumOfElements>shuffleLengthThreshold){

			RadixShuffleSort(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, mShiftRightAmounts[stage_index + 1], stage_index + 1, inputArrayIsDestination);
			AddEliminator(begin + curOffset, begin + curOffset + curNumOfElements, eliminator);
			
        }
        else{
			_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold>(begin + curOffset, followBegin + curOffset, curNumOfElements, straightshiftRightAmount);
			InsertionSortWithMoveAndEliminator(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, eliminator, !inputArrayIsDestination);
        }

    }
	else {
		if (curNumOfElements > InsertThreshold){

			introsort_detail::IntrosortRec(begin + curOffset, begin + curOffset + curNumOfElements, introsort_detail::IntrosortDepth(begin + curOffset, begin + curOffset + curNumOfElements),
				std::less<index_type>(), internal::DualSwapper<RandomIt, FollowIt>(begin + curOffset, followBegin + curOffset));
		}
		InsertionSortWithMove(begin + curOffset, followBegin + curOffset, begin2 + curOffset, followBegin2 + curOffset, curNumOfElements, !inputArrayIsDestination);
		
	}
	
	

}

template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
template<class RandomIt, class FollowIt>
void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
    RadixShuffleFind(RandomIt begin, FollowIt followBegin,size_t length)
{


    

    auto curIt=begin+1;
    auto cur_divisor=mFastDivisors[0]; //current divisors
    auto cur_num_divisor=mDivisors[0];
	index_type curValue = ((unsignedType)*begin) / cur_divisor; //get current value that defines to block to find
    index_type nextValue=(curValue+1)*cur_num_divisor; //get the next value that defines the end of current block
	//
    size_t curOffset=0; //offset and number of elements in current block
    size_t curNumOfElements=1;
    index_type eliminator;
    unsigned int straightshiftRightAmount;
    unsignedType straightEliminatorMult, nextEliminatorMult;
    size_t shuffleLengthThreshold;
    getNextShifts(0,straightshiftRightAmount,straightEliminatorMult,nextEliminatorMult,shuffleLengthThreshold );
	size_t est_length = length / mMaxDims[0];
	int num_runs = estimateNumberOfRunsShuffle(est_length, 1);
	size_t t_length = est_length;
	//std::cout << " cur_num_divisor " << cur_num_divisor << std::endl;
	int num_straight_runs = estimateNumberOfRuns<Threshold, InsertThreshold, PowerOfTwoRadix>(t_length, mTotalMax/mMaxDims[0]);


    //std::cout << " nextMaxSize " << nextMaxSize << " newStraightMaxSize " << newStraightMaxSize << std::endl;
	//std::cout << " nextEliminatorMult " << nextEliminatorMult << " straightEliminatorMult " << straightEliminatorMult << std::endl;
    //std::cout << "next_bit_size " << next_bit_size << " number_shifts " << number_shifts << " straightshiftRightAmount " << straightshiftRightAmount << " shuffleLengthThreshold " << shuffleLengthThreshold << std::endl;

    for(;curIt<begin+length;++curIt){

        if(*curIt>=nextValue){
			//std::cout << *curIt << " numelemenst " << curNumOfElements << std::endl;
            if ( curNumOfElements >= Threshold ){
				

				if (num_runs<num_straight_runs)
                    eliminator=curValue*nextEliminatorMult ;
                else
                    eliminator=curValue*straightEliminatorMult;
				//std::cout << "eliminator " << eliminator << std::endl;
                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                    *it-=eliminator;

				if (num_runs<num_straight_runs){
					//std::cout << "Met the rp threshold" << std::endl;
					
					mResize(2*curNumOfElements);
                    RadixShuffleSort(begin+curOffset,followBegin+curOffset,mScratch1.get(),mScratch2.get(),curNumOfElements,mShiftRightAmounts[1],1,false);	
					////std::cout << "Shufflesorted" << std::endl;
					AddEliminator(begin + curOffset, begin + curOffset + curNumOfElements, eliminator);
					////std::cout << "Eliminator" << std::endl;
                }
                else{
					//std::cout << "Met the straight radix sort threshold" << std::endl;
					_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold>(begin + curOffset, followBegin + curOffset, curNumOfElements, straightshiftRightAmount);
					InsertionSortImprovedEliminator(begin + curOffset, begin + curOffset + curNumOfElements, followBegin + curOffset, eliminator);
                }
				

            }
            else {
				//std::cout << "Have to intro/insert" << std::endl;
				if (curNumOfElements > InsertThreshold){
					//std::cout << "Met the intro thershold" << std::endl;
					introsort_detail::IntrosortRec(begin + curOffset, begin + curOffset + curNumOfElements, introsort_detail::IntrosortDepth(begin + curOffset, begin + curOffset + curNumOfElements),
						std::less<index_type>(), internal::DualSwapper<RandomIt, FollowIt>(begin + curOffset, followBegin + curOffset));
				}
				if (curNumOfElements>1)

					InsertionSortImproved(begin + curOffset, begin + curOffset + curNumOfElements, followBegin + curOffset, std::less<index_type>());
            }

            curOffset=curIt-begin;
            curNumOfElements=1;

            curValue=((unsignedType)*curIt)/cur_divisor;
            nextValue=(curValue+1)*cur_num_divisor;
			
        }
        else{
            ++curNumOfElements;
        }


    }

    //perform last section
    if ( curNumOfElements >= Threshold )
    {
		
		if (num_runs<num_straight_runs)
            eliminator=curValue*nextEliminatorMult ;
        else
            eliminator=curValue*straightEliminatorMult;
		//std::cout << "eliminator " << eliminator << std::endl;
        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
            *it-=eliminator;

		if (num_runs<num_straight_runs){
			//std::cout << "Radix shuff curNumOfElements " << curNumOfElements << " " << curOffset << std::endl;
			
			//_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold>(begin + curOffset, followBegin + curOffset, curNumOfElements, straightshiftRightAmount);
			//InsertionSortImprovedEliminator(begin + curOffset, begin + curOffset + curNumOfElements, followBegin + curOffset, eliminator);
			
			mResize(2 * curNumOfElements);
            RadixShuffleSort(begin+curOffset,followBegin+curOffset,mScratch1.get(),mScratch2.get(),curNumOfElements,mShiftRightAmounts[1],1,false);	
			AddEliminator(begin + curOffset, begin + curOffset + curNumOfElements, eliminator);
        }
        else{
			//std::cout << "Straight Radix sort curNumOfElements " << curNumOfElements << " " << curOffset << std::endl;
			_RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold>(begin + curOffset, followBegin + curOffset, curNumOfElements, straightshiftRightAmount);
			InsertionSortImprovedEliminator(begin + curOffset, begin + curOffset + curNumOfElements, followBegin + curOffset, eliminator);
			
        }
		

    }
    else{
		if (curNumOfElements > InsertThreshold)    {


			introsort_detail::IntrosortRec(begin + curOffset, begin + curOffset + curNumOfElements, introsort_detail::IntrosortDepth(begin + curOffset, begin + curOffset + curNumOfElements),
				std::less<index_type>(), internal::DualSwapper<RandomIt, FollowIt>(begin + curOffset, followBegin + curOffset));
		}
		if (curNumOfElements>1)
			InsertionSortImproved(begin + curOffset, begin + curOffset + curNumOfElements, followBegin + curOffset, std::less<index_type>());
		

    }
	
	


}

template<class index_type, class data_type, unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix, unsigned int Threshold, unsigned int InsertThreshold>
template<class RandomIt1, class FollowIt1, class RandomIt2, class FollowIt2>
inline void RadixShuffle<index_type, data_type, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold>::
performRadixExtractDigitRotate(RandomIt1  begin1, FollowIt1  followBegin1, RandomIt2  begin2, FollowIt2  followBegin2, const size_t length, const unsigned int cur_shiftRightAmount,const int stage_index){

	
	/*if (isFinal)
		std::cout << "Enetered Final Extract Digit" <<std::endl;*/
	libmia_constexpr unsigned int numberOfBins = PowerOfTwoRadix;
	const unsigned int total_shift_amount = cur_shiftRightAmount + mDivisorsShiftAmount[stage_index]; //the total shift amount - standard radix shift plus the shift approximation to the divisor


	//constexpr unsigned int BitMask =PowerOfTwoRadix-1;
	for (auto &i : mCountBuffer)
		i = 0;
	index_type min_index(std::numeric_limits<index_type>::max());
	index_type cur_Digit;
	auto end = begin1 + length;
	for (auto it = begin1; it < end; ++it){ // Scan the array and count the number of times each value in the current bitMask appears
		cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it), total_shift_amount);
		mCountBuffer[cur_Digit+ 1]++;
		min_index = std::min(*it, min_index);
	}
	cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)min_index), total_shift_amount); //extract the approximate min digit, which we'll use to start our scan of when to rotate the extracted digits

	for (auto idx = 1; idx<numberOfBins + 1; ++idx) // Turn the counts into starting offsets into the array
		mCountBuffer[idx] += mCountBuffer[idx - 1];

	auto followIt = followBegin1;
	auto mCountBuffer2 = mCountBuffer;
	
	
	
	
	

	//if (cur_shiftRightAmount){
		for (auto it = begin1 ; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
			cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it), total_shift_amount);
			*(begin2 + mCountBuffer2[cur_Digit]) = *it;
			*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		}

		if (!cur_shiftRightAmount){
			auto check_divisor = 1;
			if (stage_index+1<mDivisors.size())
				check_divisor=mDivisors[stage_index+1];
			//std::cout << " check divisor" << check_divisor << std::endl;
			for (auto idx = 0; idx < numberOfBins ; ++idx){				
				auto followIt = followBegin2 + mCountBuffer[idx];
				for (auto it = begin2 + mCountBuffer[idx]; it < begin2 + mCountBuffer[idx + 1]-1; ++it, ++followIt){
					if (*it>*(it+1) && *it - *(it + 1) >= check_divisor){
						//if (idx == 1){
							//std::cout << "Rotating " << *(it + 1)<< " because " << *it << " is bigger, at position " << it + 1 - begin2 << std::endl;
							//std::cout << *it - *(it + 1) << " " << mDivisors[stage_index + 1] << " " << check_divisor << std::endl;
						//}
						std::rotate(begin2 + mCountBuffer[idx], it + 1, begin2 + mCountBuffer[idx + 1]);
						std::rotate(followBegin2 + mCountBuffer[idx], followIt + 1, followBegin2 + mCountBuffer[idx + 1]);
						
						break;
					}

				}
			}

		}
	//}
	//else{
	//	unsignedType trigger = ((unsignedType)min_index / mFastDivisors[stage_index]);// >> cur_shiftRightAmount;//this is the true min digit 
	//	//std::cout << "trigger quotient of " << min_index << " divided by " << mDivisors[stage_index] << " is " << trigger << std::endl;
	//	unsignedType true_inc = mDivisors[stage_index];// << cur_shiftRightAmount;
	//	trigger = (trigger + 1)*true_inc;
	//	unsignedType start = min_index >> total_shift_amount;
	//	//std::cout << "start quotient of " << min_index << " divided by " << mDivisorsShiftAmount[stage_index] << " is " << start << std::endl;
	//	start = (start + 1) << total_shift_amount ;//-1;
	//	unsignedType inc = unsignedType(1) << total_shift_amount;
	//	
	//	std::array<index_type, numberOfBins> bufferRotated;		
	//	for (int i = 0; i < cur_Digit; ++i){
	//		bufferRotated[i] = 0;
	//	}		
	//	//std::cout << "Trigger " << trigger;
	//	for (int i = cur_Digit; i < numberOfBins; ++i){
	//		//std::cout << "bin " << i << " start " << start << std::endl;
	//		if (start > trigger){

	//			bufferRotated[i] = trigger;
	//			trigger += true_inc;
	//			//std::cout << "Triggered new value " << trigger << std::endl;
	//		}
	//		else
	//			bufferRotated[i] = 0;

	//		start += inc;
	//	}

	//	for (auto it = begin1 ; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
	//		cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it), total_shift_amount);
	//		if (bufferRotated[cur_Digit] && *it < bufferRotated[cur_Digit]){ //if we need to rotate shift right
	//			/*if (cur_Digit == 1){
	//				std::cout << "About to rotate 3 because " << *it << " < " << bufferRotated[cur_Digit] << std::endl;
	//				for (auto i = mCountBuffer[cur_Digit]; i < mCountBuffer2[cur_Digit]; ++i){
	//					std::cout << *(begin2 + i) << std::endl;
	//				}
	//			}*/
	//			//*(begin2 + mCountBuffer2[cur_Digit] - 1) = -1;
	//			//shift existing entries in the buffer right
	//			std::copy(std::reverse_iterator<RandomIt2>(begin2 + mCountBuffer2[cur_Digit]), std::reverse_iterator<RandomIt2>(begin2 + mCountBuffer[cur_Digit]), std::reverse_iterator<RandomIt2>(begin2 + mCountBuffer[cur_Digit+1]));
	//			std::copy(std::reverse_iterator<FollowIt2>(followBegin2 + mCountBuffer2[cur_Digit]), std::reverse_iterator<FollowIt2>(followBegin2 + mCountBuffer[cur_Digit]), std::reverse_iterator<FollowIt2>(followBegin2 + mCountBuffer[cur_Digit + 1]));
	//			
	//			/*if (cur_Digit == 1){
	//				std::cout << "Rotated 3 because " << *it << " < " << bufferRotated[cur_Digit] << std::endl;
	//				for (auto i = mCountBuffer[cur_Digit + 1] - (mCountBuffer2[cur_Digit] - mCountBuffer[cur_Digit]); i < mCountBuffer[cur_Digit + 1]; ++i){
	//					std::cout << *(begin2 + i) << std::endl;
	//				}
	//			}*/
	//			
	//			//reset start of cur buffer
	//			mCountBuffer2[cur_Digit] = mCountBuffer[cur_Digit];
	//			bufferRotated[cur_Digit] = 0;
	//		}
	//		
	//		*(begin2 + mCountBuffer2[cur_Digit]) = *it;
	//		*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;
	//		//if (cur_Digit == 3){
	//		//	//std::cout << "Added to 3" << std::endl;
	//		//	/*for (auto i = mCountBuffer[cur_Digit]; i < mCountBuffer2[cur_Digit]; ++i){
	//		//		std::cout << *(begin2 + i) << std::endl;
	//		//	}*/
	//		//}
	//		
	//	}
	//}

	//std::cout << "********************Final******************** " << std::endl;
	//for (auto i = 0; i < mCountBuffer.size() - 1; ++i){
	/*for (auto i = 1; i < 2; ++i){
		for (auto it = begin2 + mCountBuffer[i]; it < begin2 + mCountBuffer[i+1]; ++it){
			std::cout << *(it) << std::endl;
			if (extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / mFastDivisors[stage_index], 0) < extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*(it - 1)) / mFastDivisors[stage_index], 0)){
				std::cout << "Issue in i " << i << " digit " << extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / mFastDivisors[stage_index], 0) << std::endl;
			}
		}
	}*/

	


}


template<class index_type, class data_type, unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix, unsigned int Threshold, unsigned int InsertThreshold>
template<class RandomIt1, class FollowIt1, class RandomIt2, class FollowIt2>
void RadixShuffle<index_type, data_type, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold>::
RadixShuffleSort(RandomIt1  begin1, FollowIt1  followBegin1, RandomIt2  begin2, FollowIt2  followBegin2, size_t length,
unsigned int cur_shiftRightAmount, int stage_index, bool inputArrayIsDestination){




	using namespace std::chrono;
	typedef std::chrono::duration<float> float_seconds;
	high_resolution_clock::time_point t1, t2;
    //std::cout << "Entered radix sort - length " << length << " cur_shiftRightAmount " << cur_shiftRightAmount << " shiftRightAmounts[stage_index] " << mShiftRightAmounts[stage_index] << std::endl;
    bool stageDone=false;
	//t1 = high_resolution_clock::now();
	//std::cout << "cur_shiftRightAmount " << cur_shiftRightAmount << std::endl;
	performRadixExtractDigitRotate(begin1, followBegin1, begin2, followBegin2, length, cur_shiftRightAmount,stage_index);
	/*std::cout << mCountBuffer[1] << std::endl;
	for (auto it = begin2; it < begin2 + mCountBuffer[1]; ++it)
		std::cout << *it << std::endl;
	std::cout << "mDivisors" << mDivisors[stage_index + 1] << std::endl;*/
	//std::copy(begin2, begin2 + length, begin1);
	//return;
    if(cur_shiftRightAmount==0){
        if(stage_index+2>=mFastDivisors.size()){//if we have no more stages to perform, so just return, with a copy if needed
			//std::cout << "STage done " << std::endl;
			if (!inputArrayIsDestination){
				//std::cout << "Doing copy " << std::endl;
				std::copy(begin2, begin2 + length, begin1);
				std::copy(followBegin2, followBegin2 + length, followBegin1);
			}
			return;
        }
        else{ //otherwise, we need to proceed to the next stage, which will be a find
			stage_index++;
            stageDone=true;
        }
    }
    else{
        cur_shiftRightAmount -= Log2ofPowerOfTwoRadix; //should always be a multiple of Log2ofPowerOfTwoRadix
    }
	
	

    auto count=mCountBuffer; //free up buffer for for use in lower recursions
	inputArrayIsDestination = !inputArrayIsDestination;
	size_t offset;
	size_t numOfElements;
	//t2 = high_resolution_clock::now();
	//std::cout << "Total Radix Move Time " << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
	//t1 = high_resolution_clock::now();
	if (stageDone){
		
		
		for (int i = 0; i < count.size() - 1; ++i)
		{
			offset = count[i];
			numOfElements = count[i + 1] - count[i];
			
				


			
			/*introsort_detail::IntrosortRec(begin2 + offset, begin2 + offset + numOfElements, introsort_detail::IntrosortDepth(begin2 + offset, begin2 + offset + numOfElements),
				std::less<index_type>(), internal::DualSwapper<RandomIt2, FollowIt2>(begin2 + offset, followBegin2 + offset));
			InsertionSortWithMove(begin2 + offset, followBegin2 + offset, begin1 + offset, followBegin1 + offset, numOfElements, !inputArrayIsDestination);*/
			if (!numOfElements){
				continue;
			}
			else if (numOfElements >= Threshold){
				//std::cout << "Finding " << numOfElements << " " << offset << std::endl;
				/*if (i == 0){
					auto it2 = begin1;
					for (auto it = begin2 + offset; it < begin2 + offset + numOfElements; ++it, ++it2)
						std::cout << *it << " " << *it2 << std::endl;
				}*/
				//return;
				RadixShuffleFind(begin2 + offset, followBegin2 + offset, begin1 + offset, followBegin1 + offset, numOfElements, stage_index, inputArrayIsDestination); //radix shufflefind will perform final insertion sort cleanup

			}
			else{
				if (numOfElements > InsertThreshold){

					introsort_detail::IntrosortRec(begin2 + offset, begin2 + offset + numOfElements, introsort_detail::IntrosortDepth(begin2 + offset, begin2 + offset + numOfElements),
						std::less<index_type>(), internal::DualSwapper<RandomIt2, FollowIt2>(begin2 + offset, followBegin2 + offset));
				}
				//std::cout << "Insert Sorting " << numOfElements << " " << offset << std::endl;
				///if (inputArrayIsDestination){

					InsertionSortWithMove(begin2 + offset, followBegin2 + offset, begin1 + offset, followBegin1 + offset, numOfElements, !inputArrayIsDestination); //if we don't call the above, we need to perform our own cleanup
				//}
			}
			
		}
	}
	else{
		
		for (int i = 0; i < count.size() - 1; ++i)
		{
			offset = count[i];
			numOfElements = count[i + 1] - count[i];
			//std::cout << "Sorting again " << numOfElements << std::endl;
			if (!numOfElements)
				continue;
			else if (numOfElements >= Threshold){
				/*introsort_detail::IntrosortRec(begin2 + offset, begin2 + offset + numOfElements, introsort_detail::IntrosortDepth(begin2 + offset, begin2 + offset + numOfElements),
				std::less<index_type>(), internal::DualSwapper<RandomIt2, FollowIt2>(begin2 + offset, followBegin2 + offset));
				InsertionSortWithMove(begin2 + offset, followBegin2 + offset, begin1 + offset, followBegin1 + offset, numOfElements, !inputArrayIsDestination);
				*/
				//std::cout << "Sorting " << *(begin2 + offset) << " " << *(begin2 + offset + numOfElements-1) << std::endl;
				//std::cout << "Double shuffle" << std::endl;
				RadixShuffleSort(begin2 + offset, followBegin2 + offset, begin1 + offset, followBegin1 + offset, numOfElements, cur_shiftRightAmount, stage_index, inputArrayIsDestination);

			}
			else {
				//std::cout << "Introing" << std::endl;
				if (numOfElements > InsertThreshold){
					
					introsort_detail::IntrosortRec(begin2 + offset, begin2 + offset + numOfElements, introsort_detail::IntrosortDepth(begin2 + offset, begin2 + offset + numOfElements),
						std::less<index_type>(), internal::DualSwapper<RandomIt2, FollowIt2>(begin2 + offset, followBegin2 + offset));
				}
				//if (inputArrayIsDestination){
					InsertionSortWithMove(begin2 + offset, followBegin2 + offset, begin1 + offset, followBegin1 + offset, numOfElements, !inputArrayIsDestination);

				//}
				
			}
		}
	}
	//t2 = high_resolution_clock::now();
	//std::cout << "Next stage Time " << duration_cast<float_seconds>(t2 - t1).count() << std::endl;
}



template<class index_type, class data_type, unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix, unsigned int Threshold, unsigned int InsertThreshold>
template<class RandomIt1, class FollowIt1, class RandomIt2, class FollowIt2>
inline void RadixShuffle<index_type, data_type, PowerOfTwoRadix, Log2ofPowerOfTwoRadix, Threshold, InsertThreshold>::
performRadixExtractDigit(RandomIt1  begin1, FollowIt1  followBegin1, RandomIt2  begin2, FollowIt2  followBegin2, const size_t length, const unsigned int cur_shiftRightAmount, const int stage_index){

	using namespace libdivide;
	//std::cout << "Enetered Extract Digit" <<std::endl;
	libmia_constexpr unsigned int numberOfBins = PowerOfTwoRadix;
	auto const fast_divisor = mFastDivisors[stage_index];


	//constexpr unsigned int BitMask =PowerOfTwoRadix-1;
	for (auto &i : mCountBuffer)
		i = 0;

	auto end = begin1 + length;
	switch (fast_divisor.get_algorithm()) {
	case 0:
		for (auto it = begin1; it<end; ++it) // Scan the array and count the number of times each value in the current bitMask appears
			mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<0>(fast_divisor), cur_shiftRightAmount) + 1]++;
		break;
	case 1:
		for (auto it = begin1; it<end; ++it) // Scan the array and count the number of times each value in the current bitMask appears
			mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<1>(fast_divisor), cur_shiftRightAmount) + 1]++;
		break;
	case 2:
		for (auto it = begin1; it<end; ++it) // Scan the array and count the number of times each value in the current bitMask appears
			mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<2>(fast_divisor), cur_shiftRightAmount) + 1]++;
		break;
	case 3:
		for (auto it = begin1; it<end; ++it) // Scan the array and count the number of times each value in the current bitMask appears
			mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<3>(fast_divisor), cur_shiftRightAmount) + 1]++;
		break;
	case 4:
		for (auto it = begin1; it<end; ++it) // Scan the array and count the number of times each value in the current bitMask appears
			mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<4>(fast_divisor), cur_shiftRightAmount) + 1]++;
		break;
	}





	for (auto idx = 1; idx<numberOfBins + 1; ++idx) // Turn the counts into starting offsets into the array
		mCountBuffer[idx] += mCountBuffer[idx - 1];

	auto followIt = followBegin1;
	auto  mCountBuffer2 = mCountBuffer;


	switch (fast_divisor.get_algorithm()) {
	case 0:
	{
		auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*begin1) / unswitch<0>(fast_divisor), cur_shiftRightAmount);
		*(begin2 + mCountBuffer2[cur_Digit]) = *begin1;
		*begin1 = -length;
		*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		for (auto it = begin1 + 1; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
			cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<0>(fast_divisor), cur_shiftRightAmount);
			*(begin2 + mCountBuffer2[cur_Digit]) = *it;
			*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		}
		break;
	}

	case 1:
	{

		auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*begin1) / unswitch<1>(fast_divisor), cur_shiftRightAmount);
		*(begin2 + mCountBuffer2[cur_Digit]) = *begin1;
		*begin1 = -length;
		*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		for (auto it = begin1 + 1; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
			cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<1>(fast_divisor), cur_shiftRightAmount);
			*(begin2 + mCountBuffer2[cur_Digit]) = *it;
			*it = -1;
			*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		}
		break;
	}
	case 2:
	{

		auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*begin1) / unswitch<2>(fast_divisor), cur_shiftRightAmount);
		*(begin2 + mCountBuffer2[cur_Digit]) = *begin1;
		*begin1 = -length;
		*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;
		for (auto it = begin1 + 1; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
			cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<2>(fast_divisor), cur_shiftRightAmount);
			*(begin2 + mCountBuffer2[cur_Digit]) = *it;
			*it = -1;
			*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		}
		break;
	}
	case 3:
	{
		auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*begin1) / unswitch<3>(fast_divisor), cur_shiftRightAmount);
		*(begin2 + mCountBuffer2[cur_Digit]) = *begin1;
		*begin1 = -length;
		*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;
		for (auto it = begin1 + 1; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
			cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<3>(fast_divisor), cur_shiftRightAmount);
			*(begin2 + mCountBuffer2[cur_Digit]) = *it;
			*it = -1;
			*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		}
		break;
	}
	case 4:
	{
		auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*begin1) / unswitch<4>(fast_divisor), cur_shiftRightAmount);
		*(begin2 + mCountBuffer2[cur_Digit]) = *begin1;
		*begin1 = -length;
		*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;
		for (auto it = begin1 + 1; it < end; ++it){ //perform radix sort on the current bitMask, using the starting offsets
			cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<4>(fast_divisor), cur_shiftRightAmount);
			*(begin2 + mCountBuffer2[cur_Digit]) = *it;
			*(followBegin2 + mCountBuffer2[cur_Digit]++) = *followIt++;

		}
		break;
	}
	}





}



} // namespace internal

}// namespace LibMIA

#endif // LIBMIA_RADIX_H
