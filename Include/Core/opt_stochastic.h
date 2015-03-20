//
//  opt_stochastic.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_opt_stochastic_h
#define Metaphor_opt_stochastic_h

namespace im
{
    // Stochastic gradient function minimizer
    
    // Derive a class from this solver and implement the function to evaluate f(x) and df(x)/dx
    // Call init with a pre-allocated vector to hold the value of x
    // Remember to initialize it to your starting point guess
    // Call step() in a while loop until it returns true, i.e while(!step());
    // You can print the progress by adding print statements in the while loop that examine fx(), delta_fx(), delta_x()
    // Get the result by calling state()
    // If you want to evaluate over a validation set then you can do this inside the while loop at infrequent intervals
    // or else in eval_end_step()
    
    enum StochasticUpdateMode
    {
        StochasticUpdateModeSimple,       // x' = x - rate * gradient(x)
        StochasticUpdateModeMomentum,     // v' = momentum * v - rate * gradient(x), x' = x + v'
        StochasticUpdateModeAccelerated,  // v' = momentum * v - rate * gradient(x + momentum * v), x' = x + v'
    };
    
    struct StochasticMinParams
    {
        StochasticMinParams()
        {
            update_mode = StochasticUpdateModeSimple;
            blend_mode = false; // if true, then rate is pre-multiplied by (1-momentum) in update formulas with momentum
            fx_eval_interval = 1; // how many step calls to count between fx evaluations, i.e. how many minibatches between cost evaluations
            delta_fx_min = 1e-8; // if delta_fx falls below this threshold then we quit
            iterations_max = 1000000; // number step calls before we quit
            rate_init = 0.01; // the initial value of the learning rate (needs setting based on your application!)
            update_rate = false; // whether to auto adjust the learning rate
            rate_multiplier = 0.1; // ratio of new to old learning rate at each update
            rate_interval = 1000; // number of steps between each update
            momentum_init = 0.9; // the initial value of the momentum if used (ignored if update_momentum is true)
            update_momentum = false; // whether to auto adjust the momentum
            momentum_max = 0.999; // maximum value for momentum
            momentum_rate = 250; // momentum = min(1-0.5/(t/momentum_rate+1),momentum_max), where t is the iteration
            compute_delta_x = true; // whether to evaluate delta_x
        }
        
        StochasticUpdateMode update_mode;
        int fx_eval_interval;
        int iterations_max;
        double delta_fx_min;
        double rate_init;
        double momentum_init;
        bool update_rate;
        double rate_multiplier;
        int rate_interval;
        bool update_momentum;
        double momentum_max;
        double momentum_rate;
        bool compute_delta_x;
        bool blend_mode;
    };
    
    template <typename TT>
    class StochasticMin
    {
    public:
        StochasticMin(StochasticMinParams const &p) : m_params(p) {}
        StochasticMin() {}
        virtual ~StochasticMin() {}
        
        // If you dont want to use the default parameters, set them before calling init()
        void set_parameters(StochasticMinParams const &p) { m_params = p; }
        
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
        
        // get and set the learning rate
        TT rate() const {return m_rate; }
        void rate(TT r) { m_rate = r; }
        
        // get and set the momentum
        TT momentum() const { return m_momentum; }
        void momentum(TT m) { m_momentum = m; }
        
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
        Vec<TT> m_vgradient;
        Vec<TT> m_vintegrator;
        StochasticMinParams m_params;
        int m_iterations;
        TT m_rate;
        TT m_momentum;
    };
}


#endif
