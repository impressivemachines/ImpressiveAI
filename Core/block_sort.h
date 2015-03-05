//
//  block_sort.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/25/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_block_sort_h
#define Metaphor_block_sort_h

namespace im
{    
    // Sorting
    // Valid types are:
    //   uint8_t
    //   int16_t
    //   int
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>
    
    // Efficient in place sort for standard vectors
    // Complex comaprison is done on real part, or magnitude if abs used
    
    template <typename TT> void core_vector_sort(std::vector<TT> &v, SortDirection direction = SortDirectionAscending);
    
    // Sort the absolute magnitude
    template <typename TT> void core_vector_sort_abs(std::vector<TT> &v, SortDirection direction = SortDirectionAscending);
    
    
    // Sorting is defined using raster order from top left to bottom right
    // Complex comaprison is done on real part, or magnitude if abs used
    
    // This can sort matrix or row or column vector data
    template <typename TT> void core_block_sort(VecView<TT> vdst, SortDirection direction = SortDirectionAscending);
    template <typename TT> void core_block_sort(MtxView<TT> mdst, SortDirection direction = SortDirectionAscending);
    
    // Sort the absolute magnitude
    template <typename TT> void core_block_sort_abs(VecView<TT> vdst, SortDirection direction = SortDirectionAscending);
    template <typename TT> void core_block_sort_abs(MtxView<TT> mdst, SortDirection direction = SortDirectionAscending);
    
    
    // Sort to MtxLoc list - the actual elements are not changed
    // Complex comaprison is done on real part, or magnitude if abs used
    
    // The output index list is resized if necessary
    template <typename TT> void core_block_sort(std::vector<int> &index_list, VecView<TT> const &msrc, SortDirection direction = SortDirectionAscending);
    template <typename TT> void core_block_sort(std::vector<MtxLoc> &index_list, MtxView<TT> const &msrc, SortDirection direction = SortDirectionAscending);
    
    // Sort the absolute magnitude
    template <typename TT> void core_block_sort_abs(std::vector<int> &index_list, VecView<TT> const &msrc, SortDirection direction = SortDirectionAscending);
    template <typename TT> void core_block_sort_abs(std::vector<MtxLoc> &index_list, MtxView<TT> const &msrc, SortDirection direction = SortDirectionAscending);
    
    

    // Find the location of an extreme value
    // Complex comaprison is done on real part, or magnitude if abs used
    // Valid types are:
    //   float
    //   double
    //   std::complex<float>
    //   std::complex<double>

    template <typename TT> int core_block_location_of_max(VecView<TT> const &msrc);
    template <typename TT> MtxLoc core_block_location_of_max(MtxView<TT> const &msrc);
    
    template <typename TT> int core_block_location_of_max_abs(VecView<TT> const &msrc);
    template <typename TT> MtxLoc core_block_location_of_max_abs(MtxView<TT> const &msrc);
    
    template <typename TT> int core_block_location_of_min(VecView<TT> const &msrc);
    template <typename TT> MtxLoc core_block_location_of_min(MtxView<TT> const &msrc);
    
    template <typename TT> int core_block_location_of_min_abs(VecView<TT> const &msrc);
    template <typename TT> MtxLoc core_block_location_of_min_abs(MtxView<TT> const &msrc);
}


#endif
