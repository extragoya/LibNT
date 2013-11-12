// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.

//This file includes helper functions to create commonly used MIAs, typically defined by some implicit function

#ifndef LIBMIAHELPERS_H_INCLUDED
#define LIBMIAHELPERS_H_INCLUDED

#include <type_traits>
#include <iostream>
#include <algorithm>


#include "LibMIAUtil.h"
#include "ImplicitMIA.h"
#include "IndexUtil.h"
#include "SparseMIA.h"


//\defgroup
namespace LibMIA
{




//!Create an implicit MIA of just ones
template<class data_type,class... Dims>
auto create_ones(Dims...dims)->ImplicitMIA<data_type,sizeof...(Dims)>{

    typedef ImplicitMIA<data_type,sizeof...(Dims)> retType;
    static_assert(internal::check_mia_constructor<retType,Dims...>::type::value,"Number of dimensions must be same as <order> and each given range must be convertible to <index_type>, i.e., integer types.");
    typedef typename internal::index_type<retType>::type index_type;

    retType ret({dims...});

    auto _function=[](index_type idx){
        return data_type(1);
    };
    ret.get_function()=_function;
    return ret;
}


//!Creates a SparseMIA based on a data func that can return the data at any valid index in mia. index_func must specify how to get the next dok_index from a given dok_index
template<class data_type, size_t _order,class DataFunc, class IndexFunc, class index_type>
void make_sparse_mia_explict_from_func(SparseMIA<data_type,_order> & mia, const DataFunc & data_func, const IndexFunc& index_func,index_type _start,index_type _end){

    index_type curIdx=_start;
    while(curIdx<_end){
        mia.push_back(data_func(curIdx),curIdx); //push_back current data, index pair
        curIdx=index_func(curIdx); //get next index
    }

}

//!creates a O(h) forward difference operator. At the last row, backward difference is performed. Second index is meant to be the inner product index
template<class data_type>
SparseMIA<data_type,2> create_forward_diff(size_t m){

    typedef SparseMIA<data_type,2> retType;
    retType ret(m,m);

    ret.reserve(2*m); //diagonal and off-diagonal
    typedef typename internal::index_type<retType>::type index_type;

    auto data_func=[m,&ret](index_type idx){
        auto indices=ret.ind2sub(idx); //get full indices
        if(indices[0]==m-1){ //if we are on the last row
            if(indices[1]-indices[0]==0) //if we are the last index, return 1,
                return 1;
            else if(indices[1]-indices[0]==-1) //if we are on the off-diagonal, return -1
                return -1;
            else
                return 0;
        }
        else{ //all other rows
            if(indices[1]-indices[0] >1 || indices[1]-indices[0]<0) //if we are not on the diagonal or off-diagonal
                return 0;
            else if (indices[1]==indices[0])
                return -1;
            else
                return 1;
        }
    };

    auto index_func=[m,&ret](index_type idx){
        auto indices=ret.ind2sub(idx);
        if(indices[0]<m-2){ //any before second-last row
            if(indices[1]==indices[0])
                indices[1]++;
            else
                indices[0]++;
        }
        else if(indices[0]==m-2){ //second last row needs to set up last row for backwards difference
            if(indices[1]==indices[0])
                indices[1]++;
            else{
                indices[0]++;
                indices[1]--;
            }
        }
        else{ //last row (do backward difference)
            if(indices[1]==indices[0]-1)
                indices[1]++;
            else
                indices[0]++; //make sure our linear index is greater than dimensionality
        }
        return ret.sub2ind(indices);
    };
    make_sparse_mia_explict_from_func(ret,data_func,index_func,index_type(0),ret.dimensionality());
    assert(ret.size()==2*m);
    return ret;
}

//!creates a O(h) backward difference operator. At the first row, forward difference is performed. Second index is meant to be the inner product index
template<class data_type>
SparseMIA<data_type,2> create_backward_diff(size_t m){

    typedef SparseMIA<data_type,2> retType;
    retType ret(m,m);

    ret.reserve(2*m); //diagonal and off-diagonal
    typedef typename internal::index_type<retType>::type index_type;

    auto data_func=[m,&ret](index_type idx){
        auto indices=ret.ind2sub(idx); //get full indices
        if(indices[0]==0){ //if we are on the first row
            if(indices[1]-indices[0]==0) //if we are the first index, return -1,
                return -1;
            else if(indices[1]-indices[0]==1) //if we are on the off-diagonal, return 1
                return 1;
            else
                return 0;
        }
        else{ //all other rows
            if(indices[0]-indices[1] >1 || indices[0]-indices[1]<0) //if we are not on the diagonal or off-diagonal
                return 0;
            else if (indices[1]==indices[0])
                return 1;
            else
                return -1;
        }
    };

    auto index_func=[m,&ret](index_type idx){
        auto indices=ret.ind2sub(idx);
        if(indices[0]>0){ //any row after the first
            if(indices[1]==indices[0]-1)
                indices[1]++;
            else
                indices[0]++;
        }
        else{ //first row
            if(indices[1]==indices[0])
                indices[1]++;
            else{
                indices[0]++;
                indices[1]--;
            }
        }
        return ret.sub2ind(indices);
    };
    make_sparse_mia_explict_from_func(ret,data_func,index_func,index_type(0),ret.dimensionality());
    assert(ret.size()==2*m);
    return ret;
}

template<class data_type,size_t _order>
SparseMIA<data_type,_order> create_delta(size_t m){



    typedef SparseMIA<data_type,_order> retType;
    typedef typename internal::index_type<retType>::type index_type;
    std::array<index_type,_order> retDims;
    std::fill(retDims.begin(),retDims.end(),m);
    retType ret(retDims);

    ret.reserve(m); //diagonal and off-diagonal


    auto data_func=[](index_type idx){
        return 1;
    };

    auto index_func=[&ret](index_type idx){
        auto indices=ret.ind2sub(idx);
        for (index_type& i : indices )
        {
            i++; // increments the value in the indices
        }
        return ret.sub2ind(indices);
    };
    make_sparse_mia_explict_from_func(ret,data_func,index_func,index_type(0),ret.dimensionality());
    assert(ret.size()==m);
    return ret;


}

//!Emulates a vector MIA by creating an order+1 MIA where the new index indexes the original MIA
template<class data_type,size_t order>
SparseMIA<data_type,order+1> concat_sparse_arrays( SparseMIA<data_type,order>& a,  SparseMIA<data_type,order>& b){
    if(a.dims()!=b.dims())
        throw MIAParameterException("Indices of two SparseMIAs must be the same to concatenate them together");
    typedef SparseMIA<data_type,order+1> retType;
    typedef typename internal::index_type<retType>::type index_type;
    std::array<index_type,order+1> retDims;
    a.reset_linIdx_sequence();
    a.sort();
    b.reset_linIdx_sequence();
    b.sort();

    for(size_t ctr=0;ctr<order;++ctr)
        retDims[ctr]=a.dim(ctr);
    retDims[order]=2; //last dim is 2, to index a or b

    retType ret(retDims);
    ret.reserve(a.size()+b.size());
    auto data_it=a.data_begin();
    for(auto it=a.index_begin();it<a.index_end();++it,++data_it){
        ret.push_back(*data_it,*it);
    }
    data_it=b.data_begin();
    for(auto it=b.index_begin();it<b.index_end();++it,++data_it){
        ret.push_back(*data_it,*it+a.dimensionality());
    }
    return ret;



}




}//namespace LibMIA
#endif // LIBMIAHELPERS_H_INCLUDED
