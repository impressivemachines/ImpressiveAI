//
//  matrix.cpp
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/29/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"

// Resizes the matrix and loses any data
template <typename TT>
void im::Mtx<TT>::resize(int nrows, int ncols)
{
    IM_CHECK_LOWER_BOUNDS(nrows, 0);
    IM_CHECK_LOWER_BOUNDS(ncols, 0);
    
    if(m_mem && m_mem.unique())
        m_mem->resize(nrows*ncols);
    else
        m_mem.reset(new std::vector<TT>(nrows*ncols));
    
    m_view.wrap(nrows, ncols, ncols, 1, m_mem->data());
}

// Stops sharing the matrix by creating own private copy and guarentees packed row major layout
template <typename TT>
void im::Mtx<TT>::stop_sharing()
{
    if(!m_mem || !m_mem.unique() || (cols()>1 && col_stride()!=1) || (rows()>1 && row_stride()!=cols()))
    {
        // Ensure simple array
        std::vector<TT> *pv = new std::vector<TT>(count());
        if(m_view.is_valid())
            core_block_copy(MtxView<TT>(rows(), cols(), cols(), 1, pv->data()), m_view);
        m_mem.reset(pv);
        m_view.wrap(rows(), cols(), cols(), 1, m_mem->data());
    }
}

// Add one or more rows at the bottom of the matrix
template <typename TT>
void im::Mtx<TT>::push_back(MtxView<TT> const &mav)
{
    stop_sharing(); // in case mtx does not own its memory (makes a local copy)
    
    if(rows()==0 || cols()==0)
    {
        // Empty matrix
        resize(mav.rows(), mav.cols());
        core_block_copy(m_view, mav);
        return;
    }
    
    m_mem->resize(count() + cols() * mav.rows()); // non destructive
    
    int last = rows();
    m_view.wrap(last + mav.rows(), cols(), cols(), 1, m_mem->data());
    
    int ncols = std::min(cols(), mav.cols());
    int i,j;
    for(i=0; i<mav.rows(); i++)
    {
        for(j=0; j<ncols; j++)
            at(last + i, j) = mav(i,j);
        for(; j<cols(); j++)
            at(last + i, j) = (TT)0.0;
    }
}

// Add an new row at the bottom of the matrix and sets its element(s) to the given value
template <typename TT>
void im::Mtx<TT>::push_back(TT const &f)
{
    stop_sharing(); // in case mtx does not own memory
    
    if(rows()==0 || cols()==0)
    {
        // Empty matrix
        resize(1, 1);
        at(0,0) = f;
        return;
    }
    
    m_mem->resize(count() + cols());
    
    int last = rows();
    m_view.wrap(last + 1, cols(), cols(), 1, m_mem->data());
    
    for(int i=0; i<cols(); i++)
        at(last, i) = f;
}

// Remove a row from the bottom of the matrix
template <typename TT>
void im::Mtx<TT>::pop_back()
{
    if(rows()>0 && cols()>0)
    {
        stop_sharing(); // in case mtx does not own memory
        
        m_mem->resize(count() - cols());
        
        m_view.wrap(rows()-1, cols(), cols(), 1, m_mem->data());
    }
}

template <typename TT>
void im::Mtx<TT>::random_uniform(Rand &rnd, TT const &low, TT const &high)
{
    for(int row=0; row<rows(); row++)
        for(int col=0; col<cols(); col++)
            at(row,col) = (TT)rnd.uniform_real(core_real(low), core_real(high));
}

template <typename TT>
void im::Mtx<TT>::random_gaussian(Rand &rnd, TT const &mean, TT const &stddev)
{
    for(int row=0; row<rows(); row++)
        for(int col=0; col<cols(); col++)
            at(row,col) = stddev * (TT)rnd.gauss() + mean;
}

template <typename TT>
TT im::Mtx<TT>::sample_bilinear(float row, float col) const
{
    if(row<0 || col<0 || row>rows()-1 || col>cols()-1)
        return (TT)0;
    
    int r0 = (int)std::floor(row);
    int c0 = (int)std::floor(col);
    
    int r1 = std::min(r0+1, rows()-1);
    int c1 = std::min(c0+1, cols()-1);
    
    float ar = row - r0;
    float ac = col - c0;
    float arac = ar * ac;
    
    TT v00 = at(r0,c0);
    TT v01 = at(r0,c1);
    TT v10 = at(r1,c0);
    TT v11 = at(r1,c1);
    
    return v00 + TT(ac) * (v01 - v00) + TT(ar) * (v10 - v00) + TT(arac) * (v00 - v01 - v10 + v11);
}

static void core_priv_clamp_ranges(int &a, int &b, int &c, int lower, int upper)
{
    if(a<lower)
        a = lower;
    else if(a>upper)
        a = upper;
    if(b<lower)
        b = lower;
    else if(b>upper)
        b = upper;
    if(c<lower)
        c = lower;
    else if(c>upper)
        c = upper;
}

template <typename TT>
TT im::Mtx<TT>::sample_bicubic(float row, float col) const
{
    int nrows = rows();
    int ncols = cols();
    
    if(row<0 || col<0 || row>nrows-1 || col>ncols-1)
        return (TT)0;
    
    int r1 = (int)std::floor(row);
    int c1 = (int)std::floor(col);
    int r0 = r1 - 1;
    int c0 = c1 - 1;
    int r2 = r1 + 1;
    int c2 = c1 + 1;
    int r3 = r1 + 2;
    int c3 = c1 + 2;
    
    if(r1<1 || r1>nrows-3)
        core_priv_clamp_ranges(r0, r2, r3, 0, nrows-1);
    if(c1<1 || c1>ncols-3)
        core_priv_clamp_ranges(c0, c2, c3, 0, ncols-1);
    
    float ar = row - r1;
    float ac = col - c1;
    float mar = (1-ar);
    float mac = (1-ac);
    
    TT krow0 = (TT)((1/6.0f) * (mar * mar - 1.0f) * mar);
    TT krow1 = (TT)(0.5f * (ar * mar + 2.0f) * mar);
    TT krow2 = (TT)(0.5f * (mar * ar + 2.0f) * ar);
    TT krow3 = (TT)((1/6.0f) * (ar * ar - 1.0f) * ar);
    
    TT kcol0 = (TT)((1/6.0f) * (mac * mac - 1.0f) * mac);
    TT kcol1 = (TT)(0.5f * (ac * mac + 2.0f) * mac);
    TT kcol2 = (TT)(0.5f * (mac * ac + 2.0f) * ac);
    TT kcol3 = (TT)((1/6.0f) * (ac * ac - 1.0f) * ac);
    
    const TT *prow0 = ptr(r0,0);
    const TT *prow1 = ptr(r1,0);
    const TT *prow2 = ptr(r2,0);
    const TT *prow3 = ptr(r3,0);
    
    TT v0 = kcol0 * prow0[c0] + kcol1 * prow0[c1] + kcol2 * prow0[c2] + kcol3 * prow0[c3];
    TT v1 = kcol0 * prow1[c0] + kcol1 * prow1[c1] + kcol2 * prow1[c2] + kcol3 * prow1[c3];
    TT v2 = kcol0 * prow2[c0] + kcol1 * prow2[c1] + kcol2 * prow2[c2] + kcol3 * prow2[c3];
    TT v3 = kcol0 * prow3[c0] + kcol1 * prow3[c1] + kcol2 * prow3[c2] + kcol3 * prow3[c3];
    
    return krow0 * v0 + krow1 * v1 + krow2 * v2 + krow3 * v3;
}



template <typename TT>
im::Mtx<TT> im::Mtx<TT>::inverse() const
{
    IM_CHECK_VALID(m_view);
    
    Mtx<TT> matRtn(rows(), cols());
    
    if(rows()!=cols() || rows()==0 || cols()==0)
    {
        matRtn = (TT)0;
        return matRtn;
    }
    
    int size = rows();
    
    if(size==1)
    {
        TT A00 = at(0,0);
        
        if(std::abs(A00) < TypeProperties<TT>::epsilon())
        {
            // Singular
            matRtn = (TT)0;
            return matRtn;
        }
        
        matRtn(0,0) = (TT)1/A00;
        
        return matRtn;
    }
    
    if(size==2)
    {
        TT det = at(0,0) * at(1,1) - at(0,1) * at(1,0);
        
        if(std::abs(det) < TypeProperties<TT>::epsilon())
        {
            // Singular
            matRtn = (TT)0;
            return matRtn;
        }
        
        det = (TT)1 / det;
        
        matRtn(0,0) = det * at(1,1);
        matRtn(0,1) = -det * at(0,1);
        matRtn(1,0) = -det * at(1,0);
        matRtn(1,1) = det * at(0,0);
        
        return matRtn;
    }
    
    // Gaussian elimination
    matRtn.set_identity();
    Mtx<TT> matSrc(size,size);
    matSrc.copy_from(m_view);
    
    for(int col = 0; col < size; col++)
    {
        int pivot_row = -1;
        double pivot_max = 0.0;
        
        for(int row = col; row < size; row++)
        {
            double f = std::abs(matSrc(row,col));
            if(f>pivot_max)
            {
                pivot_max = f;
                pivot_row = row;
            }
        }
        
        if(pivot_max < TypeProperties<TT>::epsilon())
        {
            // Singular
            matRtn = (TT)0;
            return matRtn;
        }
        
        if(pivot_row != col)
        {
            // Exchange rows
            core_block_exchange(matSrc.view().row(pivot_row), matSrc.view().row(col));
            core_block_exchange(matRtn.view().row(pivot_row), matRtn.view().row(col));
        }
        
        TT pivot_scale = ((TT)1.0) / matSrc(col,col);
        
        core_block_scale(matSrc.view().row(col), matSrc.view().row(col), pivot_scale);
        core_block_scale(matRtn.view().row(col), matRtn.view().row(col), pivot_scale);
        
        for(int row = 0; row < size; row++)
            if(row != col)
            {
                TT f = -matSrc(row,col);
                for(int i=0; i<size; i++)
                {
                    matSrc(row, i) += f * matSrc(col, i);
                    matRtn(row, i) += f * matRtn(col, i);
                }
            }
    }
    
    return matRtn;
}

template <typename TT>
TT im::Mtx<TT>::det() const
{
    MatrixDecompLU<TT> lu(m_view);
    return lu.det();
}

template <typename TT>
TT im::Mtx<TT>::log_det() const
{
    MatrixDecompLU<TT> lu(m_view);
    return lu.log_det();
}

template <typename TT>
TT im::Mtx<TT>::sign_det() const
{
    MatrixDecompLU<TT> lu(m_view);
    return lu.sign_det();
}

enum TokenType
{
    TokenTypeOpen, TokenTypeClose, TokenTypeNumber, TokenTypeRow, TokenTypeError
};

struct Token
{
    Token(TokenType t) : type(t), value(im::Cd(0)) {}
    TokenType type;
    im::Cd value;
    
    void get(float &v) { v = value.real(); }
    void get(double &v) { v = value.real(); }
    void get(im::Cf &v) { v = im::Cf(value.real(), value.imag()); }
    void get(im::Cd &v) { v = value; }
};

static Token core_priv_get_token(char const *&ptext)
{
    // skip space
    while(*ptext && isspace(*ptext))
        ptext++;
    
    if(*ptext=='[')
    {
        ptext++;
        return Token(TokenTypeOpen);
    }
    if(*ptext=='\0')
        return Token(TokenTypeClose);
    if(*ptext==']')
    {
        ptext++;
        return Token(TokenTypeClose);
    }
    if(*ptext==';')
    {
        ptext++;
        return Token(TokenTypeRow);
    }
    
    char *ptmp;
    double valre = std::strtod(ptext, &ptmp);
    if(ptmp==ptext)
        return Token(TokenTypeError);
    
    ptext = ptmp;
    
    // skip space
    while(*ptext && isspace(*ptext))
        ptext++;
    
    Token tok(TokenTypeNumber);
    
    if(*ptext=='i')
    {
        ptext++;
        tok.value = im::Cd(0, valre);
        return tok;
    }
    
    tok.value = im::Cd(valre, 0);
    
    char const *psave = ptext;
    
    if(*ptext=='-' || *ptext=='+')
    {
        // possible complex number
        int sign = 1;
        if(*ptext=='-')
            sign = -1;
        
        ptext++;
        
        double valim = std::strtod(ptext, &ptmp);
        if(ptext==ptmp)
            return Token(TokenTypeError);
        
        ptext = ptmp;
        
        // skip space
        while(*ptext && isspace(*ptext))
            ptext++;
        
        if(*ptext=='i')
        {
            ptext++;
            tok.value = im::Cd(valre, valim * sign);
        }
        else
            ptext = psave;
    }
    
    return tok;
}

// parse matlab format data
template <typename TT>
im::Mtx<TT> &im::Mtx<TT>::operator<<(char const *ptext)
{
    std::vector<TT> vecval;
    std::vector<MtxLoc> vecloc;
    
    IM_CHECK_NULL(ptext);
    
    int cols = 0;
    int r = 0;
    int c = 0;
    
    bool rowcomplete = false;
    bool started = false;
    
    while(true)
    {
        Token tok = core_priv_get_token(ptext);
        switch(tok.type)
        {
            case TokenTypeOpen:
                if(started)
                {
                    IM_THROW_BAD_DATA;
                    resize(0, 0);
                    return *this;
                }
                break;
                
            case TokenTypeClose:
            {
                // wrap up
                if(r==0 && cols==0)
                {
                    resize(0, 0);
                    return *this;
                }
                
                resize(r+1, cols);
                set((TT)0);
                for(int i=0; i<vecloc.size(); i++)
                {
                    MtxLoc &loc = vecloc[i];
                    at(loc.row, loc.col) = vecval[i];
                }
                
                return *this;
            }
                
            case TokenTypeNumber:
            {
                if(rowcomplete)
                {
                    c = 0;
                    r++;
                    rowcomplete = false;
                }
                
                TT val;
                tok.get(val);
                vecval.push_back(val);
                vecloc.push_back(MtxLoc(r,c));
                c++;
                if(c>cols)
                    cols = c;
                break;
            }
                
            case TokenTypeRow:
                rowcomplete = true;
                break;
                
            case TokenTypeError:
            default:
                IM_THROW_BAD_DATA;
                resize(0, 0);
                return *this;
        }
        
        started = true;
    }
}

// Returns index list for values matching comparison
template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_equal_to(TT const &val)
{
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(at(i,j)==val)
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_greater_than(TT const &val)
{
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(core_real(at(i,j))>core_real(val))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_less_than(TT const &val)
{
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(core_real(at(i,j))<core_real(val))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_greater_than_abs(TT const &val)
{
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(std::abs(at(i,j)) > std::abs(val))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_less_than_abs(TT const &val)
{
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(std::abs(at(i,j)) < std::abs(val))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

// Returns index list for values matching comparison against elements of equal size matrix

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_equal_to(Mtx const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_MATRIX_SIZES_MATCH(m_view, src);
    
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(at(i,j)==src(i,j))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_greater_than(Mtx const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_MATRIX_SIZES_MATCH(m_view, src);
    
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(core_real(at(i,j)) > core_real(src(i,j)))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_less_than(Mtx const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_MATRIX_SIZES_MATCH(m_view, src);
    
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(core_real(at(i,j)) < core_real(src(i,j)))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_greater_than_abs(Mtx const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_MATRIX_SIZES_MATCH(m_view, src);
    
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(std::abs(at(i,j)) > std::abs(src(i,j)))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

template <typename TT>
std::vector<im::MtxLoc> im::Mtx<TT>::select_less_than_abs(Mtx const &src)
{
    IM_CHECK_VALID(src);
    IM_CHECK_MATRIX_SIZES_MATCH(m_view, src);
    
    std::vector<MtxLoc> index;
    for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
            if(std::abs(at(i,j)) < std::abs(src(i,j)))
            {
                MtxLoc loc(i,j);
                index.push_back(loc);
            }
    
    return index;
}

// Instantiate classes
template class im::Mtx<float>;
template class im::Mtx<double>;
template class im::Mtx<im::Cf>;
template class im::Mtx<im::Cd>;
