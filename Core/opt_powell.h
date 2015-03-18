//
//  opt_powell.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_opt_powell_h
#define Metaphor_opt_powell_h

namespace im
{
    // Powell's direction set minimizer is a gradient-free method of finding a minimum
    // of a function of N variables.
    
    // Derive a class from this solver and implement the function to evaluate f(x)
    // Call init with a pre-allocated vector to hold the value of x
    // Remember to initialize it to your starting point guess
    // Call step() in a while loop until it returns true, i.e while(!step());
    // You can print the progress by adding print statements in the while loop that examine fx(), delta_fx(), delta_x()
    // Get the result by calling state()
    
    struct PowellMinParams
    {
        PowellMinParams()
        {
            bracket_max = 1.0; // bracketing operation in line minimization uses the range [0, bracket_max]
            termination_ratio = 1e-7; // ratio of delta_fx to fx to cause termination
            iterations_max = 1000000; // number step calls before we quit
            line_min_eps = 0; // if >0 set the error tolerance for line min termination
        }
        
        double bracket_max;
        double termination_ratio;
        int iterations_max;
        double line_min_eps;
    };

    // the actual minimizer
    template <typename TT>
    class PowellMin
    {
    public:
        PowellMin(PowellMinParams const &p) : m_params(p) {}
        PowellMin() {}
        virtual ~PowellMin() {}
        
        // If you dont want to use the default parameters, set them before calling init()
        void set_parameters(PowellMinParams const &p) { m_params = p; }
        
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
        // The only essential function is eval_fx.
        
        // Called after init is called
        virtual void eval_init() {}
        
        // Called at the start of each step.
        virtual void eval_start_step() {}
        
        // Compute the function f(X), given vector X.
        virtual TT eval_fx(Vec<TT> const &vx) = 0;
        
        // not used by powell
        virtual void eval_dfx(Vec<TT> &vdfx, Vec<TT> const &vx) { vdfx = (TT)0; }
        
        // Called at the end of each step. Return true if you need to early exit.
        virtual bool eval_end_step() { return false; }
        
    protected:
        bool m_early_exit;
        bool m_startup;
        TT m_fx;
        TT m_delta_fx;
        TT m_delta_x;
        Vec<TT> m_vstate;
        Vec<TT> m_vstatesave;
        Mtx<TT> m_mdirset;
        Vec<TT> m_vdir;
        PowellMinParams m_params;
        VectorLineMin<TT, PowellMin<TT> > m_linemin;
        int m_iterations;
    };
    
}

#endif
