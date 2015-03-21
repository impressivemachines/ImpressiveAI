//
//  opt_vstepgrad.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/18/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_opt_vstepgrad_h
#define Metaphor_opt_vstepgrad_h

namespace im
{
    // THIS PARTICULAR CODE IS NOT YET COMPLETE MARCH21/2015
    // Simon
    
    // Variable stepsize gradient descent function minimizer
    
    // Derive a class from this solver and implement the function to evaluate f(x) and df(x)/dx
    // Call init with a pre-allocated vector to hold the value of x
    // Remember to initialize it to your starting point guess
    // Call step() in a while loop until it returns true, i.e while(!step());
    // You can print the progress by adding print statements in the while loop that examine fx(), delta_fx(), delta_x()
    // Get the result by calling state()
    
    struct VariableStepGradientMinParams
    {
        VariableStepGradientMinParams()
        {
            iterations_max = 1000000; // number step calls before we quit
            termination_ratio = 1e-7; // ratio of delta_fx to fx to cause termination
            stepsize_init = 0.001; // initial value for the stepsize
            turning_threshold = 0.9; // minimum dot product to accept a gradient update
            stepsize_change_on_accept = 1.5; // step multiplier if gradient turning slowly (must be > 1)
            stepsize_change_on_reject = 0.66; // step multplier if gradient changing too quickly (must be < 1)
            momentum = 0.9; // gradient IIR smoothing filter coefficient
        }
    
        int iterations_max;
        double termination_ratio;
        double stepsize_init;
        double turning_threshold;
        double stepsize_change_on_accept;
        double stepsize_change_on_reject;
        double momentum;
        
    };
    
    template <typename TT>
    class VariableStepGradientMin
    {
    public:
        VariableStepGradientMin(VariableStepGradientMinParams const &p) : m_params(p) {}
        VariableStepGradientMin() {}
        virtual ~VariableStepGradientMin() {}
        
        // If you dont want to use the default parameters, set them before calling init()
        void set_parameters(VariableStepGradientMinParams const &p) { m_params = p; }
        
        // Initalize.
        // You can wrap external memory or pass in a freshly allocated vector
        void init(Vec<TT> vstate);
        
        // Take one minimization step. Call this in a loop until it returns true
        // Then get the result by calling state()
        bool step();
        
        // Returns true if the solver stopped before convergence.
        bool early_exit() const { return m_early_exit; }
        
        // Current value of the function
        TT fx() const { return m_fx; }
        
        // Get the change in function value during the last step
        TT delta_fx() const { return m_delta_fx; }
        
        // Get the distance moved during the last step (note step sizes may rise and fall)
        TT delta_x() const { return m_delta_x; }
        
        // Get/set the current state vector x
        Vec<TT> state() { return m_vstate; }
        
        // Dimenstionality of the state vector x
        int dims() const { return m_vstate.rows(); }
        
        // Returns the number of calls to step()
        int iteration_count() const { return m_iterations; }
        
        // Functions which are to be over-ridden by derived class.
        
        // Called after init is called
        virtual void eval_init() {}
        
        // Called at the start of each step.
        virtual void eval_start_step() {}
        
        // Compute the function f(X), given vector X.
        // this is where you would evaluate the cost fuction over the whole dataset
        // this information is only used to test for termination, so you can return zero without
        // affecting the algorithm if you want to terminate it manually or after a fixed number of steps
        // in which case you should set delta_fx_min to -1
        // you can also use it to log the cost function over the validation set
        virtual TT eval_fx(Vec<TT> const &vx) = 0;
        
        // Compute the gradient vector of f(X), about X
        // this is where you would evaluate the mean gradient for a minibatch
        virtual void eval_dfx(Vec<TT> &vdfx, Vec<TT> const &vx) = 0;
        
        // Called at the end of each step. Return true if you need to early exit.
        virtual bool eval_end_step() { return false; }
        
    protected:
        bool m_early_exit;
        bool m_startup;
        TT m_fx;
        TT m_delta_fx;
        TT m_delta_x;
        Vec<TT> m_vstate;
        Vec<TT> m_vtrialstate;
        Vec<TT> m_vgrad;
        Vec<TT> m_vtrialgrad;
        VariableStepGradientMinParams m_params;
        int m_iterations;
        TT m_stepsize;
        //TT m_gradmag;
        //TT m_trialgradmag;
    };

}

#endif
