//
//  optimize.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_optimize_h
#define Metaphor_optimize_h

namespace im
{
    template <typename TT> class MinBase;
    
    // Base class for function minimization
    // Derive your function class from this one
    template <typename TT>
    class FuncEval
    {
    public:
        // Called by optimizer to initialize the vector x
        virtual void eval_init(VecView<TT> vx) {  vx = (TT)0; }
        
        // Called by optimizer at the start of a new step
        virtual void eval_start_step(MinBase<TT> &opt) {}
        
        // Evaluate the function at x
        // e.g.
        // return core_block_reduce_sum_squares(vx);
        virtual TT eval_fx(VecView<TT> const &vx) { return (TT)0; }
        
        // Evaluate the derivative at x
        // e.g.
        // for(int i=0; i<vx.rows(); i++)
        //    vdf(i) = 2*vx(i);
        virtual void eval_dfx(VecView<TT> vdf, VecView<TT> const &vx) { IM_THROW_NO_IMPL; }
        
        // Evaluate the hessian at x
        // The hessian is a symmetric matrix. Only the lower triangular portion should be set. Use the
        // TRI_LO(N,i,j) indexing macro to set the elements.
        // e.g.
        // for(int i=0; i<vx.rows(); i++)
        //    for(int j=0; j<=i; j++)
        //       vhess(TRI_LO(vx.rows(), i, j)) = ((i==j) ? (TT)2 : (TT)0);
        virtual void eval_ddfx(VecView<TT> vhess, VecView<TT> const &vx) { IM_THROW_NO_IMPL; }
        
        // Called by the optimizer at the end of each step
        // If it returns false then the optimizer will exit early
        virtual bool eval_end_step(MinBase<TT> &opt) { return true; }
    };

    // Virtual base class for function minimizing
    template <typename TT>
    class MinBase
    {
    public:
        // The dimensionality of the state vector is indicated by the number of rows of vstate.
        // The memory for the state vector can be allocated by the caller or by the optimizer.
        // Set vstate = VecView<TT>(n,0,NULL) to let the optimizer do the allocation.
        // Before it returns, this function initializes the state vector by calling pfe->eval_init()
        virtual void init(FuncEval<TT> *pfe, VecView<TT> vstate) = 0;
        
        // To optimize, call step() repeatedly in a loop until it returns false.
        // This means that either it is finished, or else the FuncEval class terminated the optimization early.
        // Step first calls pfe->eval_start_step(), then does optimization which involves calling one or more of the
        // eval functions, and finally calls pfe->eval_end_step(), and then it returns.
        virtual bool step() = 0;
        
        // Determine if the optimizer actually completed.
        virtual bool is_finished() = 0;
        
        // Report the distance that was moved over the last step
        virtual TT delta_fx() = 0;
        virtual TT delta_x() = 0;
        
        // Get a reference to the state vector (x).
        virtual VecView<TT> &state() = 0;
    };
    
    

}

#endif
