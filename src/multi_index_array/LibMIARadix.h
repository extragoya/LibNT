
#ifndef LIBMIA_RADIX_H
#define LIBMIA_RADIX_H
#include "LibMIAAlgorithm.h"
namespace LibMIA
{







namespace internal
{

// Copyright(c), Victor J. Duvanenko, 2009

// Listing 3
template<unsigned int BitMask, class _Type1>
inline typename std::make_unsigned<_Type1>::type extractDigit( _Type1 a, unsigned int shiftRightAmount )
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
template<unsigned long PowerOfTwoRadix ,unsigned long Log2ofPowerOfTwoRadix, long Threshold, class ForwardIt, class FollowIt >
inline void _RadixSort_Unsigned_PowerOf2Radix_L1( ForwardIt begin, FollowIt followBegin, size_t length, unsigned long shiftRightAmount )
{


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
            _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix, Threshold >( begin+ startOfBin[ i ],followBegin+startOfBin[i], numberOfElements , shiftRightAmount );
        else if ( numberOfElements >= 50 )
            introsort_detail::IntrosortRec(begin+ startOfBin[ i ],begin+ startOfBin[ i ]+numberOfElements,introsort_detail::IntrosortDepth(begin+ startOfBin[ i ],begin+ startOfBin[ i ]+numberOfElements),
                    std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin+ startOfBin[ i ],followBegin+ startOfBin[ i ]));

    }

}



template< class ForwardIt,class FollowIt>
inline void RadixSortInPlace_PowerOf2Radix_Unsigned(ForwardIt  begin, FollowIt  followBegin , size_t length, size_t max_size)
{

    typedef typename std::iterator_traits<ForwardIt>::value_type index_type;

    const long Threshold                      =  3000;
    const unsigned long PowerOfTwoRadix       = 2048;    // Radix - must be a power of 2
    const unsigned long Log2ofPowerOfTwoRadix =   11;    // log( 2048 ) = 11
    auto current_bit_size=(unsigned long)std::ceil(std::log2(max_size));

    auto number_shifts=(unsigned long)std::ceil((double)current_bit_size/Log2ofPowerOfTwoRadix);

    unsigned long shiftRightAmount;
    if(current_bit_size>=Log2ofPowerOfTwoRadix)
        shiftRightAmount= (number_shifts-1)*Log2ofPowerOfTwoRadix;
    else
        shiftRightAmount=0;


    if ( length >= Threshold )
        _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix, Threshold >( begin, followBegin, length ,  shiftRightAmount );
    else
        introsort_detail::IntrosortRec(begin,begin+length,introsort_detail::IntrosortDepth(begin,begin+length),
                    std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin,followBegin));

    //clean up the last
    //InsertionSort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<ForwardIt,FollowIt>(begin,followBegin));
    InsertionSortImproved(begin,begin+length,followBegin,std::less<index_type>());
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
inline int estimateNumberOfRuns(size_t length,size_t max_size){

    int numberOfRuns=1;
    int radixRuns=static_cast<int>(std::ceil(std::log2(max_size)/std::log2(Radix)));
    size_t first_run=max_size/static_cast<size_t>(std::pow(Radix,radixRuns-1));
    //std::cout << "First Run " << first_run << " radixRuns " << radixRuns << std::endl;
    length/=first_run;
    radixRuns--;
    while(length>Threshold && radixRuns>0){
        numberOfRuns++;
        radixRuns--;
        length/=Radix;
    }
    if(radixRuns>0 && length >IntroThreshold)
        numberOfRuns++;
    return numberOfRuns;
}

template<unsigned int IntroThreshold, unsigned int Radix>
inline size_t getLengthThreshold(size_t max_size){


	int radixRuns = static_cast<int>(std::ceil(std::log2(max_size) / std::log2(Radix)));
    size_t first_run=max_size/static_cast<size_t>(std::pow(Radix,radixRuns-1));
    size_t thresholdLength=50;
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
    std::vector<index_type> mScratch1;
    std::vector<data_type> mScratch2;
    std::array<size_t,PowerOfTwoRadix+1> mCountBuffer;
    std::array<size_t,PowerOfTwoRadix+1> mCountBuffer2;
    const size_t mTotalMax;
    unsigned int mTotalBitSize;
    unsigned int mTotalRightAmount;
    const bool mFirstSortOrFind;

    typedef typename std::vector<index_type>::iterator ScratchIt1;
    typedef typename std::vector<data_type>::iterator ScratchIt2;
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
    :mMaxDims(_MaxDims),mDivisors(_Divisors),mTotalMax(_TotalMax),mFirstSortOrFind(_SortOrFind)
    {
        assert(mMaxDims.size()==mDivisors.size());
        auto number_of_stages=mMaxDims.size();
        mShiftRightAmounts.resize(mMaxDims.size());
        mFastDivisors.resize(mMaxDims.size());
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

        }
    }

    void mResize(size_t length){
        if(mScratch1.size()<length){
            mScratch1.resize(length);
            mScratch2.resize(length);
        }
    }

    //!Perform the permutation
    template<class RandomIt, class FollowIt>
    void permute(RandomIt begin, FollowIt  followBegin,size_t length)
    {
        if(length<Threshold)
            Introsort(begin,begin+length,std::less<index_type>(),internal::DualSwapper<RandomIt,FollowIt>(begin,followBegin));


        if(mFirstSortOrFind==false)  //if we are currently finding and not sorting
        {


            RadixShuffleFind(begin,followBegin,length,0);


        }
        else
        {
            //std::cout << "max size " << mMaxDims.front() << " length " << length << " estimate of runs " << estimateNumberOfRuns<Threshold,InsertThreshold,PowerOfTwoRadix>(length,mMaxDims.front()) << std::endl;
            //std::cout << "total_max_size " << mTotalMax << " estimate of runs " << estimateNumberOfRuns<Threshold,InsertThreshold,PowerOfTwoRadix>(length,mTotalMax) << std::endl;

            if (estimateNumberOfRuns<Threshold,InsertThreshold,PowerOfTwoRadix>(length,mMaxDims.front()) < estimateNumberOfRuns<Threshold,InsertThreshold,PowerOfTwoRadix>(length,mTotalMax))
            {

                mResize(length);
                RadixShuffleSort(begin,followBegin,mScratch1.begin(),mScratch2.begin(),length,mShiftRightAmounts[0],0,false);

            }
            else
            {
                //std::cout << "Doing straight radix - not worth it " << mTotalRightAmount << std::endl;
                //RadixSortInPlace_PowerOf2Radix_Unsigned( begin,  followBegin,length,total_max_size );
                _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold>( begin,followBegin, length, mTotalRightAmount );
                //RadixSortStraight<PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold>(  begin,  followBegin, scratch1.begin(),scratch2.begin(),length,total_max_size);
            }

        }
        InsertionSortImproved(begin,begin+length,followBegin,std::less<index_type>());
    }
private:
    //only should run as the first recursion
    template<class RandomIt, class FollowIt>
    void RadixShuffleFind(RandomIt begin, FollowIt followBegin,size_t length,int stage_index);

    //only should run as later recurions
    template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
    void RadixShuffleFind(RandomIt begin, FollowIt followBegin,RandomIt2 begin2, FollowIt2 followBegin2,size_t length,int stage_index,bool inputArrayIsDestination);

    template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
    void RadixShuffleSort(RandomIt begin, FollowIt followBegin,RandomIt2 begin2, FollowIt2 followBegin2,size_t length,unsigned int curShiftRightAmount,int stage_index,bool inputArrayIsDestination);

    template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
    void performRadixExtractDigit(RandomIt begin1,FollowIt followBegin1,RandomIt2 begin2,FollowIt2 followBegin2,size_t length,unsigned int cur_shiftRightAmount,int stage_index);

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


template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
template<class RandomIt, class FollowIt,class RandomIt2, class FollowIt2>
void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
    RadixShuffleFind(RandomIt begin, FollowIt followBegin,RandomIt2 begin2, FollowIt2 followBegin2,size_t length,int stage_index,bool inputArrayIsDestination)
{


    if(stage_index>=mFastDivisors.size()-1) //if this is the last stage (which is also a user error), just return, as there's no reason to do a find
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

    for(;curIt<begin+length;++curIt){

        if(*curIt>=nextValue){




            if ( curNumOfElements >= Threshold ){


                if(curNumOfElements>shuffleLengthThreshold)
                    eliminator=curValue*straightEliminatorMult;
                else
                    eliminator=curValue*nextEliminatorMult;

                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                    *it-=eliminator;

                if(curNumOfElements>shuffleLengthThreshold){

                    RadixShuffleSort(begin+curOffset,followBegin+curOffset,mScratch1.begin(),mScratch2.begin(),curNumOfElements,mShiftRightAmounts[stage_index+1],stage_index+1,inputArrayIsDestination);
                    if ( inputArrayIsDestination ){
                        for(auto it=begin2+curOffset;it< begin2+curOffset+curNumOfElements;++it)
                            *it+=eliminator;

                    }
                    else{
                        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                            *it+=eliminator;
                    }

                }
                else{

                    _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold>( begin+curOffset,followBegin+curOffset, curNumOfElements, straightshiftRightAmount);
                    if ( inputArrayIsDestination ){
                        auto it2=begin2;
                        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it,++it2)
                            *it2=*it+eliminator;
                        std::copy(followBegin+curOffset,followBegin+curOffset+curNumOfElements,followBegin2+curOffset);
                    }
                    else{
                        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                            *it+=eliminator;
                    }
                }


            }
            else if(curNumOfElements>InsertThreshold){

                introsort_detail::IntrosortRec(begin+curOffset, begin+curOffset+curNumOfElements,introsort_detail::IntrosortDepth(begin+curOffset, begin+curOffset+curNumOfElements),
                    std::less<index_type>(),internal::DualSwapper<RandomIt,FollowIt>(begin+curOffset,followBegin+curOffset));
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
        if(curNumOfElements>shuffleLengthThreshold)
            eliminator=curValue*straightEliminatorMult;
        else
            eliminator=curValue*nextEliminatorMult;

        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
            *it-=eliminator;

        if(curNumOfElements>shuffleLengthThreshold){

            RadixShuffleSort(begin+curOffset,followBegin+curOffset,mScratch1.begin(),mScratch2.begin(),curNumOfElements,mShiftRightAmounts[stage_index+1],stage_index+1,inputArrayIsDestination);
            if ( inputArrayIsDestination ){
                for(auto it=begin2+curOffset;it< begin2+curOffset+curNumOfElements;++it)
                    *it+=eliminator;

            }
            else{
                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                    *it+=eliminator;
            }

        }
        else{

            _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold>( begin+curOffset,followBegin+curOffset, curNumOfElements, straightshiftRightAmount);
            if ( inputArrayIsDestination ){
                auto it2=begin2;
                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it,++it2)
                    *it2=*it+eliminator;
                std::copy(followBegin+curOffset,followBegin+curOffset+curNumOfElements,followBegin2+curOffset);
            }
            else{
                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                    *it+=eliminator;
            }
        }

    }
    else if(curNumOfElements>InsertThreshold)    {


        introsort_detail::IntrosortRec(begin+curOffset, begin+curOffset+curNumOfElements,introsort_detail::IntrosortDepth(begin+curOffset, begin+curOffset+curNumOfElements),
                                       std::less<index_type>(),internal::DualSwapper<RandomIt,FollowIt>(begin+curOffset,followBegin+curOffset));

    }


}

template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
template<class RandomIt, class FollowIt>
void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
    RadixShuffleFind(RandomIt begin, FollowIt followBegin,size_t length,int stage_index)
{


    if(stage_index>=mFastDivisors.size()-1) //if this is the last stage (which is also a user error), just return, as there's no reason to do a find
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

    for(;curIt<begin+length;++curIt){

        if(*curIt>=nextValue){




            if ( curNumOfElements >= Threshold ){


                if(curNumOfElements>shuffleLengthThreshold)
                    eliminator=curValue*straightEliminatorMult;
                else
                    eliminator=curValue*nextEliminatorMult;

                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                    *it-=eliminator;

                if(curNumOfElements>shuffleLengthThreshold){
                    mResize(2*curNumOfElements);
                    RadixShuffleSort(begin+curOffset,followBegin+curOffset,mScratch1.begin(),mScratch2.begin(),curNumOfElements,mShiftRightAmounts[stage_index+1],stage_index+1,false);
                }
                else{

                    _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold>( begin+curOffset,followBegin+curOffset, curNumOfElements, straightshiftRightAmount);

                }
                for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
                    *it+=eliminator;


            }
            else if(curNumOfElements>InsertThreshold){

                introsort_detail::IntrosortRec(begin+curOffset, begin+curOffset+curNumOfElements,introsort_detail::IntrosortDepth(begin+curOffset, begin+curOffset+curNumOfElements),
                    std::less<index_type>(),internal::DualSwapper<RandomIt,FollowIt>(begin+curOffset,followBegin+curOffset));
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
        if(curNumOfElements>shuffleLengthThreshold)
            eliminator=curValue*straightEliminatorMult;
        else
            eliminator=curValue*nextEliminatorMult;

        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
            *it-=eliminator;

        if(curNumOfElements>shuffleLengthThreshold){
            mResize(2*curNumOfElements);
            RadixShuffleSort(begin+curOffset,followBegin+curOffset,mScratch1.begin(),mScratch2.begin(),curNumOfElements,mShiftRightAmounts[stage_index+1],stage_index+1,false);
        }
        else{

            _RadixSort_Unsigned_PowerOf2Radix_L1<PowerOfTwoRadix,Log2ofPowerOfTwoRadix,Threshold>( begin+curOffset,followBegin+curOffset, curNumOfElements, straightshiftRightAmount);
        }
        for(auto it=begin+curOffset;it< begin+curOffset+curNumOfElements;++it)
            *it+=eliminator;

    }
    else if(curNumOfElements>InsertThreshold)    {


        introsort_detail::IntrosortRec(begin+curOffset, begin+curOffset+curNumOfElements,introsort_detail::IntrosortDepth(begin+curOffset, begin+curOffset+curNumOfElements),
                                       std::less<index_type>(),internal::DualSwapper<RandomIt,FollowIt>(begin+curOffset,followBegin+curOffset));

    }


}

template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
template<class RandomIt1, class FollowIt1,class RandomIt2, class FollowIt2>
inline void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
performRadixExtractDigit( RandomIt1  begin1, FollowIt1  followBegin1,RandomIt2  begin2, FollowIt2  followBegin2, size_t length,unsigned int cur_shiftRightAmount,int stage_index){

    using namespace libdivide;
    //std::cout << "Enetered Extract Digit" <<std::endl;
    libmia_constexpr unsigned int numberOfBins=PowerOfTwoRadix;
    auto fast_divisor=mFastDivisors[stage_index];


    //constexpr unsigned int BitMask =PowerOfTwoRadix-1;
    for(auto &i:mCountBuffer)
        i=0;


   switch (fast_divisor.get_algorithm()) {
        case 0:
            for ( auto it=begin1;it<begin1+length;++it ) // Scan the array and count the number of times each value in the current bitMask appears
				mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<0>(fast_divisor), cur_shiftRightAmount) + 1]++;
            break;
        case 1:
            for ( auto it=begin1;it<begin1+length;++it ) // Scan the array and count the number of times each value in the current bitMask appears
				mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<1>(fast_divisor), cur_shiftRightAmount) + 1]++;
            break;
        case 2:
            for ( auto it=begin1;it<begin1+length;++it ) // Scan the array and count the number of times each value in the current bitMask appears
				mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<2>(fast_divisor), cur_shiftRightAmount) + 1]++;
            break;
        case 3:
            for ( auto it=begin1;it<begin1+length;++it ) // Scan the array and count the number of times each value in the current bitMask appears
				mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<3>(fast_divisor), cur_shiftRightAmount) + 1]++;
            break;
        case 4:
            for ( auto it=begin1;it<begin1+length;++it ) // Scan the array and count the number of times each value in the current bitMask appears
				mCountBuffer[extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<4>(fast_divisor), cur_shiftRightAmount) + 1]++;
            break;
    }





    for(auto idx=1;idx<numberOfBins+1;++idx) // Turn the counts into starting offsets into the array
        mCountBuffer[idx]+=mCountBuffer[idx-1];

    auto followIt=followBegin1;
    auto mCountBuffer2=mCountBuffer;


    switch (fast_divisor.get_algorithm()) {
        case 0:
            for ( auto it=begin1;it<begin1+length;++it,++followIt){ //perform radix sort on the current bitMask, using the starting offsets
				auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<0>(fast_divisor), cur_shiftRightAmount);
                *(begin2+mCountBuffer2[cur_Digit])=*it;
                *(followBegin2+mCountBuffer2[cur_Digit]++)=*followIt;
            }
            break;
        case 1:
            for ( auto it=begin1;it<begin1+length;++it,++followIt){ //perform radix sort on the current bitMask, using the starting offsets
				auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<1>(fast_divisor), cur_shiftRightAmount);
                *(begin2+mCountBuffer2[cur_Digit])=*it;
                *(followBegin2+mCountBuffer2[cur_Digit]++)=*followIt;
            }
            break;
        case 2:
            for ( auto it=begin1;it<begin1+length;++it,++followIt){ //perform radix sort on the current bitMask, using the starting offsets
				auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<2>(fast_divisor), cur_shiftRightAmount);
                *(begin2+mCountBuffer2[cur_Digit])=*it;
                *(followBegin2+mCountBuffer2[cur_Digit]++)=*followIt;
            }
            break;
        case 3:
            for ( auto it=begin1;it<begin1+length;++it,++followIt){ //perform radix sort on the current bitMask, using the starting offsets
				auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<3>(fast_divisor), cur_shiftRightAmount);
                *(begin2+mCountBuffer2[cur_Digit])=*it;
                *(followBegin2+mCountBuffer2[cur_Digit]++)=*followIt;
            }
            break;
        case 4:
            for ( auto it=begin1;it<begin1+length;++it,++followIt){ //perform radix sort on the current bitMask, using the starting offsets
				auto cur_Digit = extractDigit<PowerOfTwoRadix - 1>(((unsignedType)*it) / unswitch<4>(fast_divisor), cur_shiftRightAmount);
                *(begin2+mCountBuffer2[cur_Digit])=*it;
                *(followBegin2+mCountBuffer2[cur_Digit]++)=*followIt;
            }
            break;
    }





}

template<class index_type, class data_type,unsigned int PowerOfTwoRadix, unsigned int Log2ofPowerOfTwoRadix,unsigned int Threshold,unsigned int InsertThreshold>
template<class RandomIt1, class FollowIt1,class RandomIt2, class FollowIt2>
void RadixShuffle<index_type, data_type,PowerOfTwoRadix, Log2ofPowerOfTwoRadix,Threshold,InsertThreshold>::
RadixShuffleSort( RandomIt1  begin1, FollowIt1  followBegin1,RandomIt2  begin2, FollowIt2  followBegin2, size_t length,
                    unsigned int cur_shiftRightAmount,int stage_index,bool inputArrayIsDestination){





    //std::cout << "Entered radix sort - length " << length << " cur_shiftRightAmount " << cur_shiftRightAmount << " shiftRightAmounts[stage_index] " << mShiftRightAmounts[stage_index] << std::endl;
    bool stageDone=false;



    performRadixExtractDigit(begin1,followBegin1,begin2,followBegin2,length,cur_shiftRightAmount,stage_index);

    if(cur_shiftRightAmount==0){

        if(stage_index+2>=mFastDivisors.size()){//if we have no more stages to perform, so just return, with a copy if needed

            if ( !inputArrayIsDestination ){

                std::copy(begin2,begin2+length,begin1);
                std::copy(followBegin2,followBegin2+length,followBegin1);
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
    for( int i = 0; i < count.size()-1; ++i )
    {
        size_t offset=count[i];
        size_t numOfElements = count[ i+1 ] - count[ i ];

        if ( numOfElements >= Threshold ){
            //std::cout << "Doing a radix recursion i " << i << " " << offset << " " << numOfElements << std::endl;
            if(stageDone){
                RadixShuffleFind(begin2+offset,followBegin2+offset,begin1+offset,followBegin1+offset,numOfElements,stage_index,inputArrayIsDestination);
            }
            else{
                RadixShuffleSort( begin2+offset,followBegin2+offset,begin1+offset,followBegin1+offset,numOfElements , cur_shiftRightAmount, stage_index,inputArrayIsDestination);
            }
        }
        else {
            if(numOfElements>InsertThreshold){

                introsort_detail::IntrosortRec(begin2+offset,begin2+offset+numOfElements,introsort_detail::IntrosortDepth(begin2+offset,begin2+offset+numOfElements),
                        std::less<index_type>(),internal::DualSwapper<RandomIt2,FollowIt2>(begin2+offset,followBegin2+offset));
            }
            if ( inputArrayIsDestination ){ //normally this would be !inputArrayIsDestination, but we inverted the value above
                std::copy(begin2+offset,begin2+offset+numOfElements,begin1+offset);
                std::copy(followBegin2+offset,followBegin2+offset+numOfElements,followBegin1+offset);
            }
        }
    }


}






} // namespace internal

}// namespace LibMIA

#endif // LIBMIA_RADIX_H
