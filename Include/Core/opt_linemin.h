//
//  opt_linemin.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/10/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_opt_linemin_h
#define Metaphor_opt_linemin_h

namespace im
{
    // Base class for 1d function evaluation
    // Derive your function class from this
    template <typename TT>
    class FuncEval1D
    {
    public:
        virtual TT eval_fx(TT x) { return (TT)0; }
        virtual TT eval_dfx(TT x) { return (TT)0; } // optional
    };
    
    // bracket a 1D minimum starting from the interval ab and expanding out
    // throws exception if a minimum was not found
    template <typename TT>
    void core_line_min_bracket(TT &xa, TT &xb, TT &xc, FuncEval1D<TT> *peval, TT bmin, TT bmax);
    
    // find a 1D minimum without derivatives
    // throws exception if the minimum can't be found in the given interval
    template <typename TT>
    void core_line_min(TT &xmin, TT &fxmin, FuncEval1D<TT> *peval, TT bmin, TT bmax, TT eps = (TT)0);
    
    // find a 1D minimum using derivatives
    // throws exception if the minimum can't be found in the given interval
    template <typename TT>
    void core_line_min_using_derivs(TT &xmin, TT &fxmin, FuncEval1D<TT> *peval, TT bmin, TT bmax, TT eps = (TT)0);
    
    
    // class used by minimizers
    template <typename TT, typename MINI>
    class VectorLineMin : public FuncEval1D<TT>
    {
    public:
        void init(MINI *p, int dims, TT eps, bool use_derivs)
        {
            m_p = p;
            m_vx.resize(dims);
            m_vderiv.resize(dims);
            m_eps = eps;
            m_use_derivs = use_derivs;
        }
        
        void linemin(Vec<TT> &vbase, Vec<TT> const &vdir, TT bracketmax)
        {
            m_vbase = vbase; // reference
            m_vdir = vdir; // reference
            
            TT xa, xb, xc;
            core_line_min_bracket(xa, xb, xc, this, (TT)0, bracketmax);
            if(m_use_derivs)
                core_line_min_using_derivs(m_xmin, m_fxmin, this, xa, xc, m_eps);
            else
                core_line_min(m_xmin, m_fxmin, this, xa, xc, m_eps);
            
            // copy back into vbase
            core_block_blas_axpy(vbase.view(), m_vdir.view(), m_xmin);
        }
        
        TT xmin() { return m_xmin; }
        TT fxmin() { return m_fxmin; }
        
        TT eval_fx(TT x)
        {
            m_vx.copy_from(m_vbase);
            core_block_blas_axpy(m_vx.view(), m_vdir.view(), x);
            return m_p->eval_fx(m_vx);
        }
        
        TT eval_dfx(TT x)
        {
            // eval derivative
            // note that core_line_min_using_derivs always calls eval_dfx after calling eval_fx
            // and uses the same value of x, so it is not necessary to re-calculate m_vx
            m_p->eval_dfx(m_vderiv, m_vx);
            return m_vderiv.dot_product(m_vdir); // project back onto the line direction
        }
        
    private:
        MINI *m_p;
        TT m_xmin;
        TT m_fxmin;
        Vec<TT> m_vbase;
        Vec<TT> m_vdir;
        Vec<TT> m_vx;
        Vec<TT> m_vderiv;
        bool m_use_derivs;
        TT m_eps;
    };
}

#endif
