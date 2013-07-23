// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.



#ifndef FUNCTION_UTIL_H
#define FUNCTION_UTIL_H




#include "Util.h"
#include "MIA.h"
#include "IndexUtil.h"
namespace LibMIA
{

namespace internal{
/** \addtogroup util Utilities
 *  @{
*/

template<class Derived,class otherDerived, class Op, class index_param_type>
typename MIAMergeReturnType<Derived,otherDerived>::type
perform_implicit_merge(const MIA<Derived>& a, const MIA<otherDerived>& b,const Op& op,const std::array<index_param_type,internal::order<Derived>::value>& index_order){


    typedef typename MIAMergeReturnType<Derived,otherDerived>::type retMIAType;

    typedef typename internal::function_type<retMIAType>::type function_type;
    typedef typename internal::index_type<retMIAType>::type index_type;
    function_type _function;

    if(internal::check_ascending(index_order)){ //if there is no index shuffle, the function is simpler and faster

       _function=[&a,&b,op](index_type idx){
            return op(a.atIdx(idx),a.convert(b.atIdx(idx)));
            //return this->atIdx(idx)+b.atIdx(idx);
        };
    }
    else{
        _function=[&a,&b,op,index_order](index_type idx){
            return op(a.atIdx(idx),a.convert(b.atIdx(internal::sub2ind(a.ind2sub(idx),index_order,b.dims()))));
        };
    }
    return retMIAType(_function,a.dims());

}



template<typename Derived, class idx_typeR, class idx_typeC, class idx_typeT, size_t R, size_t C, size_t T>
auto latticeCopy(const MIA<Derived> &mia, const std::array<idx_typeR,R> & row_indices, const std::array<idx_typeC,C> & column_indices,const std::array<idx_typeT,T> & tab_indices)
->DenseLattice<typename internal::data_type<Derived>::type>
{

    typedef typename internal::index_type<Derived>::type index_type;
    typedef typename internal::data_type<Derived>::type data_type;
    constexpr auto order= internal::order<Derived>::value ;
    static_assert(internal::check_index_compatibility<index_type,idx_typeR>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeC>::type::value,"Must use an array convertable to index_type");
    static_assert(internal::check_index_compatibility<index_type,idx_typeT>::type::value,"Must use an array convertable to index_type");
    static_assert(R+C+T==order,"Size of all three arrays must equal mOrder");
    //statically check number of indices match up
    size_t row_size=1, column_size=1, tab_size=1;
    std::array<index_type, R> row_dims;
    std::array<index_type, C> column_dims;
    std::array<index_type, T> tab_dims;

    //std::cout <<"Tab " << tab_indices[0] << " " << tab_indices.size() << "\n";
    //std::cout <<"Dims " << this->m_dims[0] << " " << this->m_dims.size() << "\n";

    row_size=internal::reorder_from(mia.dims(), row_indices,row_dims);
    column_size=internal::reorder_from(mia.dims(), column_indices,column_dims);
    tab_size=internal::reorder_from(mia.dims(), tab_indices,tab_dims);

    //std::cout<< "Tab dims " << tab_dims[0] << "\n";
    DenseLattice<data_type> lat(row_size, column_size, tab_size);



    index_type row_idx=0, column_idx=0,tab_idx=0;
    for (index_type k=0;k<tab_size;++k){
        tab_idx=internal::sub2ind(internal::ind2sub(k,tab_dims),tab_indices,mia.dims());
        for (index_type j=0;j<column_size;++j){
            column_idx=tab_idx+internal::sub2ind(internal::ind2sub(j,column_dims),column_indices,mia.dims());
            for (index_type i=0;i<row_size;++i){
                row_idx=column_idx+internal::sub2ind(internal::ind2sub(i,row_dims),row_indices,mia.dims());

                lat(i,j,k)=mia.atIdx(row_idx);

            }
        }


    }

    //lat.print();
    return lat;


}


template<typename MIA,typename otherMIA, typename array_type,size_t Inter,size_t L_outer,size_t R_outer>
auto implicitNoLatticeMult(const MIA &a,const otherMIA &b,const std::array<array_type,Inter>&l_inter_idx,const std::array<array_type,L_outer>&l_outer_idx,const std::array<array_type,Inter>&r_inter_idx,const std::array<array_type,R_outer>&r_outer_idx)
->typename MIANoLatticeProductReturnType<MIA,otherMIA,L_outer+R_outer+Inter>::type
{

    typedef typename MIANoLatticeProductReturnType<MIA,otherMIA,L_outer+R_outer+Inter>::type RetType;
    typedef typename internal::index_type<RetType>::type index_type;
    typedef typename internal::index_type<MIA>::type a_index_type;
    typedef typename internal::index_type<otherMIA>::type b_index_type;
    typedef typename internal::function_type<RetType>::type function_type;
    static_assert(internal::check_index_compatibility<index_type,array_type>::type::value,"Must use an array convertable to index_type");
    std::array<a_index_type, Inter> l_inter_dims;
    std::array<a_index_type, L_outer> l_outer_dims;
    std::array<b_index_type, Inter> r_inter_dims;
    std::array<b_index_type, R_outer> r_outer_dims;




    //get inter and outer dimensionality and the individual dimensions that make up that number - should default to one if any of the arrays are empty
    size_t l_inter_size=internal::reorder_from(a.dims(), l_inter_idx,l_inter_dims);
    size_t r_inter_size= internal::reorder_from(b.dims(), r_inter_idx,r_inter_dims);
    if(l_inter_size!=r_inter_size || !std::equal(l_inter_dims.begin(),l_inter_dims.end(),r_inter_dims.begin()))
        throw DimensionMismatchException("Element-wise dimensions must match during MIA multiplication");
    internal::reorder_from(a.dims(), l_outer_idx,l_outer_dims);
    internal::reorder_from(b.dims(), r_outer_idx,r_outer_dims);

    std::array<index_type,L_outer+R_outer+Inter> retDims;
    concat_arrays(l_outer_dims, r_outer_dims,l_inter_dims,retDims);
    RetType c(retDims);
    //create lambda function
    function_type _function=[&a,&b,&c,l_inter_idx,r_inter_idx,l_outer_idx,r_outer_idx](index_type _index){
        auto full_indices=internal::ind2sub(_index,c.dims());
        a_index_type l_idx(0);
        b_index_type r_idx(0);
        l_idx+=internal::sub2ind(full_indices.begin(),full_indices.begin()+L_outer,l_outer_idx,a.dims());
        l_idx+=internal::sub2ind(full_indices.begin()+L_outer+R_outer,full_indices.end(),l_inter_idx,a.dims());
        r_idx+=internal::sub2ind(full_indices.begin()+L_outer,full_indices.begin()+L_outer+R_outer,r_outer_idx,b.dims());
        r_idx+=internal::sub2ind(full_indices.begin()+L_outer+R_outer,full_indices.end(),r_inter_idx,b.dims());
//        print_array(full_indices,"full_indices");
//        print_array(l_outer_idx,"l_outer_idx");
//        std::cout << "Outer l_idx " << internal::sub2ind(full_indices.begin(),full_indices.begin()+L_outer,l_outer_idx,a.dims()) << std::endl;
//        std::cout << "L_outer " << L_outer << " R_outer " << R_outer << " Inter " << Inter << std::endl;
//        std::cout << "l_idx " << l_idx << " r_idx " << r_idx << std::endl;
        return a.atIdx(l_idx)*b.atIdx(r_idx);
    };
    c.get_function()=_function;
    return c;



}

/*! @} */
}





}


#endif // FUNCTION_UTIL_H

