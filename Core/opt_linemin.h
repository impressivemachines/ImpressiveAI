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
    void core_line_min(TT &xmin, TT &fxmin, FuncEval1D<TT> *peval, TT bmin, TT bmax);
    
    // find a 1D minimum using derivatives
    // throws exception if the minimum can't be found in the given interval
    template <typename TT>
    void core_line_min_using_derivs(TT &xmin, TT &fxmin, FuncEval1D<TT> *peval, TT bmin, TT bmax);
    
}

#endif
