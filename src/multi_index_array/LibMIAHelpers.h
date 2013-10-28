// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef LIBMIAHELPERS_H_INCLUDED
#define LIBMIAHELPERS_H_INCLUDED

#include <type_traits>
#include <iostream>
#include <algorithm>


#include "LibMIAUtil.h"
#include "ImplicitMIA.h"
#include "IndexUtil.h"



//\defgroup
namespace LibMIA
{


template<class data_type,size_t _order>
ImplicitMIA<data_type,_order> create_delta(size_t dim){



    typedef typename internal::index_type<ImplicitMIA<data_type,_order>>::type index_type;
    std::array<index_type,_order> dims;
    dims.fill(dim);
    ImplicitMIA<data_type,_order> delta(dims);

    auto _function=[dims](index_type idx){
            auto indices=internal::ind2sub(idx,dims);
            auto first=indices[0];
            size_t i;
            for(i=1;i<indices.size();++i){
                if(indices[i]!=first)
                    break;
            }
            if(i==indices.size())
                return 1;
            else
                return 0;

    };

    delta.get_function()=_function;
    return delta;


}

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
}

#endif // LIBMIAHELPERS_H_INCLUDED
