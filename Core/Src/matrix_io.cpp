//
//  matrix_io.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/29/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

// Read the file header only
void im::core_read_imp_file_header(FILE *fp, ImpInfo *pinfo, std::vector<uint32_t> &vsizertn)
{
    IM_CHECK_NULL(fp);
    IM_CHECK_NULL(pinfo);
    
    if(fread(pinfo, sizeof(ImpInfo), 1, fp)!=1)
        IM_THROW_FILE_ERROR;
    
    if(core_is_big_endian())
    {
        core_byte_reverse_in_place(pinfo->magic);
        core_byte_reverse_in_place(pinfo->version);
        core_byte_reverse_in_place(pinfo->elsize);
        core_byte_reverse_in_place(pinfo->core_type_id);
        core_byte_reverse_in_place(pinfo->color_model);
        core_byte_reverse_in_place(pinfo->dims);
    }

    if(pinfo->magic != IM_IMP_FILE_MAGIC || pinfo->version!=1 || pinfo->dims==0 || pinfo->dims > 65536)
        IM_THROW_FILE_FORMAT_ERROR;
    
    vsizertn.resize(pinfo->dims);
    
    for(int i=0; i<pinfo->dims; i++)
    {
        uint32_t tmp;
        
        if(fread(&tmp, sizeof(uint32_t), 1, fp) != 1)
            IM_THROW_FILE_ERROR;
        
        if(core_is_big_endian())
            core_byte_reverse_in_place(tmp);
        
        vsizertn[i] = tmp;
    }
}

//

// Just read the data for the matrix, assuming that it is correctly sized
template <typename TT> void im::core_read_imp_file_data(FILE *fp, VecView<TT> vv)
{
    IM_CHECK_NULL(fp);
    
    if(vv.row_stride()==1)
    {
        if(fread(vv.ptr(), sizeof(TT), vv.rows(), fp) != vv.rows())
            IM_THROW_FILE_ERROR;
    }
    else
    {
        for(int i=0; i<vv.rows(); i++)
        {
            if(fread(vv.ptr(i), sizeof(TT), 1, fp) != 1)
                IM_THROW_FILE_ERROR;
        }
    }
}

#define INST(TT) template void im::core_read_imp_file_data(FILE *fp, VecView<TT> vv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

// Just read the data for the matrix, assuming that it is correctly sized
template <typename TT> void im::core_read_imp_file_data(FILE *fp, MtxView<TT> mv)
{
    IM_CHECK_NULL(fp);
    
    int y;
    for(y=0; y<mv.rows(); y++)
    {
        if(mv.col_stride()==1)
        {
            if(fread(mv.ptr(y,0), sizeof(TT), mv.cols(), fp) != mv.cols())
                IM_THROW_FILE_ERROR;
        }
        else
        {
            int x;
            for(x=0; x<mv.cols(); x++)
            {
                if(fread(mv.ptr(y,x), sizeof(TT), 1, fp) != 1)
                    IM_THROW_FILE_ERROR;
            }
        }
    }
}

#define INST(TT) template void im::core_read_imp_file_data(FILE *fp, MtxView<TT> mv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

// Reads and checks the file header and then loads the data
// The matrix must already be the correct type and size

template <typename TT> void im::core_read_imp_file(FILE *fp, VecView<TT> const &vv)
{
    ImpInfo info;
    std::vector<uint32_t> vsizes;
    
    core_read_imp_file_header(fp, &info, vsizes);
    
    int const id = im::TypeProperties<TT>::core_type_id;
    
    if(sizeof(TT) != info.elsize || id != info.core_type_id)
        IM_THROW_FILE_FORMAT_ERROR;
    
    int c = 1;
    for(int i=0; i<vsizes.size(); i++)
        c *= vsizes[i];
    if(vv.rows()!=c)
        IM_THROW_FILE_FORMAT_ERROR;
    
    core_read_imp_file_data(fp, vv);
}

#define INST(TT) template void im::core_read_imp_file(FILE *fp, VecView<TT> const &vv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_read_imp_file(FILE *fp, MtxView<TT> const &mv)
{
    ImpInfo info;
    std::vector<uint32_t> vsizes;
    
    core_read_imp_file_header(fp, &info, vsizes);
    
    int const id = im::TypeProperties<TT>::core_type_id;
    
    if(sizeof(TT) != info.elsize || id != info.core_type_id)
        IM_THROW_FILE_FORMAT_ERROR;
    
    if(vsizes.size()==1)
    {
        if(mv.count()!=vsizes[0])
            IM_THROW_FILE_FORMAT_ERROR;
    }
    else if(vsizes.size()==2)
    {
        if(mv.rows()!=vsizes[0] || mv.cols()!=vsizes[1])
            IM_THROW_FILE_FORMAT_ERROR;
    }
    else
    {
        int c = 1;
        for(int i=1; i<vsizes.size(); i++)
            c *= vsizes[i];
        if(mv.rows()!=vsizes[0] || mv.cols()!=c)
            IM_THROW_FILE_FORMAT_ERROR;
    }
    
    core_read_imp_file_data(fp, mv);
}

#define INST(TT) template void im::core_read_imp_file(FILE *fp, MtxView<TT> const &mv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

// Reads and checks the file header and then loads the data
// The matrix must already be the correct type and size
template <typename TT> void im::core_load_imp_file(char const *pfilename, VecView<TT> vv)
{
    IM_CHECK_NULL(pfilename);
    
    FILE *fp = fopen(pfilename, "rb");
    core_read_imp_file(fp, vv);
    fclose(fp);
}

#define INST(TT) template void im::core_load_imp_file(char const *pfilename, VecView<TT> vv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_load_imp_file(char const *pfilename, MtxView<TT> mv)
{
    IM_CHECK_NULL(pfilename);
    
    FILE *fp = fopen(pfilename, "rb");
    core_read_imp_file(fp, mv);
    fclose(fp);
}

#define INST(TT) template void im::core_load_imp_file(char const *pfilename, MtxView<TT> mv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

// Write the file header only
void im::core_write_imp_file_header(FILE *fp, ImpInfo *pinfo, std::vector<uint32_t> const &vsize)
{
    IM_CHECK_NULL(fp);
    IM_CHECK_NULL(pinfo);
    IM_CHECK_ARGS(pinfo->elsize>0 && vsize.size()==pinfo->dims);
    IM_CHECK_ARGS(pinfo->dims>0 && pinfo->dims<=65536);
    
    ImpInfo infowrite = *pinfo;
    
    if(core_is_big_endian())
    {
        core_byte_reverse_in_place(infowrite.magic);
        core_byte_reverse_in_place(infowrite.version);
        core_byte_reverse_in_place(infowrite.elsize);
        core_byte_reverse_in_place(infowrite.core_type_id);
        core_byte_reverse_in_place(infowrite.color_model);
        core_byte_reverse_in_place(infowrite.dims);
    }
    
    if(fwrite(&infowrite, sizeof(ImpInfo), 1, fp)!=1)
        IM_THROW_FILE_ERROR;
    
    for(int i=0; i<pinfo->dims; i++)
    {
        uint32_t tmp = vsize[i];
        if(core_is_big_endian())
            core_byte_reverse_in_place(tmp);
        if(fwrite(&tmp, sizeof(uint32_t), 1, fp) != 1)
            IM_THROW_FILE_ERROR;
    }
}

//

// Write only the matrix data, assuming the file header was alredy written
template <typename TT> void im::core_write_imp_file_data(FILE *fp, VecView<TT> const &vv)
{
    IM_CHECK_NULL(fp);
    
    if(vv.row_stride()==1)
    {
        if(fwrite(vv.ptr(), sizeof(TT), vv.rows(), fp) != vv.rows())
            IM_THROW_FILE_ERROR;
    }
    else
    {
        int i;
        for(i=0; i<vv.rows(); i++)
        {
            if(fwrite(vv.ptr(i), sizeof(TT), 1, fp) != 1)
                IM_THROW_FILE_ERROR;
        }
    }
}

#define INST(TT) template void im::core_write_imp_file_data(FILE *fp, VecView<TT> const &vv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_write_imp_file_data(FILE *fp, MtxView<TT> const &mv)
{
    IM_CHECK_NULL(fp);
    
    int y;
    for(y=0; y<mv.rows(); y++)
    {
        if(mv.col_stride()==1)
        {
            if(fwrite(mv.ptr(y,0), sizeof(TT), mv.cols(), fp) != mv.cols())
                IM_THROW_FILE_ERROR;
        }
        else
        {
            int x;
            for(x=0; x<mv.cols(); x++)
            {
                if(fwrite(mv.ptr(y,x), sizeof(TT), 1, fp) != 1)
                    IM_THROW_FILE_ERROR;
            }
        }
    }
}

#define INST(TT) template void im::core_write_imp_file_data(FILE *fp, MtxView<TT> const &mv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

// Write both file header and matrix data
template <typename TT> void im::core_write_imp_file(FILE *fp, VecView<TT> const &vv)
{
    ImpInfo info;
    
    info.magic = IM_IMP_FILE_MAGIC;
    info.version = 1;
    info.elsize = sizeof(TT);
    info.core_type_id = im::TypeProperties<TT>::core_type_id;
    info.color_model = 0; // Undefined for matrices
    info.dims = 1;
    
    std::vector<uint32_t> vsizes;
    vsizes.resize(1);
    vsizes[0] = vv.rows();
    
    core_write_imp_file_header(fp, &info, vsizes);
    core_write_imp_file_data(fp, vv);
}

#define INST(TT) template void im::core_write_imp_file(FILE *fp, VecView<TT> const &vv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_write_imp_file(FILE *fp, MtxView<TT> const &mv)
{
    ImpInfo info;
    
    info.magic = IM_IMP_FILE_MAGIC;
    info.version = 1;
    info.elsize = sizeof(TT);
    info.core_type_id = im::TypeProperties<TT>::core_type_id;
    info.color_model = 0; // Undefined for matrices
    info.dims = 2;
    
    std::vector<uint32_t> vsizes;
    vsizes.resize(2);
    vsizes[0] = mv.rows();
    vsizes[1] = mv.cols();

    core_write_imp_file_header(fp, &info, vsizes);
    core_write_imp_file_data(fp, mv);
}

#define INST(TT) template void im::core_write_imp_file(FILE *fp, MtxView<TT> const &mv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

//

// Write both file header and matrix data
template <typename TT> void im::core_save_imp_file(char const *pfilename, VecView<TT> const &vv)
{
    IM_CHECK_NULL(pfilename);
    
    FILE *fp = fopen(pfilename, "wb");
    core_write_imp_file(fp, vv);
    fclose(fp);
}

#define INST(TT) template void im::core_save_imp_file(char const *pfilename, VecView<TT> const &vv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST

template <typename TT> void im::core_save_imp_file(char const *pfilename, MtxView<TT> const &mv)
{
    IM_CHECK_NULL(pfilename);
    
    FILE *fp = fopen(pfilename, "wb");
    core_write_imp_file(fp, mv);
    fclose(fp);
}

#define INST(TT) template void im::core_save_imp_file(char const *pfilename, MtxView<TT> const &mv)
INST(uint8_t); INST(int16_t); INST(int); INST(float); INST(double); INST(std::complex<float>); INST(std::complex<double>);
#undef INST



