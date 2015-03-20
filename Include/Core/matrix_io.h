//
//  matrix_io.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/29/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_matrix_io_h
#define Metaphor_matrix_io_h

#define IM_IMP_FILE_MAGIC             0x49504d49


namespace im
{
    // Standard header structure for imp files (for saving vectors, matrices, tensors and images)
    // After this structure the sizes of the dimensions are stored as uint32_t in major to minor order
    
    struct ImpInfo
    {
        uint32_t magic;
        uint32_t version;
        uint32_t elsize;
        uint32_t core_type_id;
        uint32_t color_model;
        uint32_t dims;
    };
    
    // Read the file header only
    void core_read_imp_file_header(FILE *fp, ImpInfo *pinfo, std::vector<uint32_t> &vsizertn);
    
    // Just read the data for the matrix, assuming that it is correctly sized
    template <typename TT> void core_read_imp_file_data(FILE *fp, VecView<TT> vv);
    template <typename TT> void core_read_imp_file_data(FILE *fp, MtxView<TT> mv);
    
    // Reads and checks the file header and then loads the data
    // The matrix must already be the correct type and size
    template <typename TT> void core_read_imp_file(FILE *fp, VecView<TT> const &vv);
    template <typename TT> void core_read_imp_file(FILE *fp, MtxView<TT> const &mv);
    
    // Reads and checks the file header and then loads the data
    // The matrix must already be the correct type and size
    template <typename TT> void core_load_imp_file(char const *pfilename, VecView<TT> vv);
    template <typename TT> void core_load_imp_file(char const *pfilename, MtxView<TT> mv);
    
    // Write the file header only
    void core_write_imp_file_header(FILE *fp, ImpInfo *pinfo, std::vector<uint32_t> const &vsize);
    
    // Write only the matrix data, assuming the file header was alredy written
    template <typename TT> void core_write_imp_file_data(FILE *fp, VecView<TT> const &vv);
    template <typename TT> void core_write_imp_file_data(FILE *fp, MtxView<TT> const &mv);
    
    // Write both file header and matrix data
    template <typename TT> void core_write_imp_file(FILE *fp, VecView<TT> const &vv);
    template <typename TT> void core_write_imp_file(FILE *fp, MtxView<TT> const &mv);
    
    // Write both file header and matrix data
    template <typename TT> void core_save_imp_file(char const *pfilename, VecView<TT> const &vv);
    template <typename TT> void core_save_imp_file(char const *pfilename, MtxView<TT> const &mv);
}

#endif
