//
//  optimize.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_optimize_h
#define Metaphor_optimize_h

// line min (powell)
// steepest descent
// stohastic
// conjugate gradient
// levenberg marquardt
// quasi newton methods
// limited memory bfgs
// graph cuts energy min

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
        void func_init(VecView<TT> vx) {  vx = (TT)0; }
        
        // Called by optimizer at the start of a new step
        void opt_start_step(MinBase<TT> &opt) {}
        
        // Evaluate the function at x
        // e.g.
        // return core_block_reduce_sum_squares(vx);
        TT func_eval(VecView<TT> const &vx) { return (TT)0; }
        
        // Evaluate the derivative at x
        // e.g.
        // for(int i=0; i<vx.rows(); i++)
        //    vdf(i) = 2*vx(i);
        void func_eval_derivative(VecView<TT> vdf, VecView<TT> const &vx) { IM_THROW_NO_IMPL; }
        
        // Evaluate the hessian at x
        // The hessian is a symmetric matrix. Only the lower triangular portion should be set. Use the
        // TRI_LO(N,i,j) indexing macro to set the elements.
        // e.g.
        // for(int i=0; i<vx.rows(); i++)
        //    for(int j=0; j<=i; j++)
        //       vhess(TRI_LO(vx.rows(), i, j)) = ((i==j) ? (TT)2 : (TT)0);
        void func_eval_hessian(VecView<TT> vhess, VecView<TT> const &vx) { IM_THROW_NO_IMPL; }
        
        // Called by the optimizer at the end of each step
        // If it returns false then the optimizer will exit early
        bool func_end_step(MinBase<TT> &opt) { return true; }
    };

    // Virtual base class for function minimizing
    template <typename TT>
    class MinBase
    {
    public:
        // The dimensionality of the state vector is indicated by the number of rows of vstate.
        // The memory for the state vector can be allocated by the caller or by the optimizer.
        // Set vstate = VecView<TT>(n,0,NULL) to let the optimizer do the allocation.
        // Before it returns, this function initializes the state vector by calling pfe->func_init()
        void init(FuncEval<TT> *pfe, VecView<TT> vstate) = 0;
        
        // To optimize, call step() repeatedly in a loop until it returns false.
        // This means that either it is finished, or else the FuncEval class terminated the optimization early.
        // Step first calls pfe->func_start_step(), then does optimization which involves calling one or more of the
        // eval functions, and finally calls pfe->func_end_step(), and then it returns.
        bool step() = 0;
        
        // Determine if the optimizer actually completed.
        bool is_finished() = 0;
        
        // Report the distance that was moved over the last step
        TT delta() = 0;
        
        // Get a reference to the state vector (x).
        VecView<TT> &state() = 0;
    };

}

#endif
