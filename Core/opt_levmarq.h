//
//  opt_levmarq.h
//  Metaphor
//
//  Created by SIMON WINDER on 3/12/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_opt_levmarq_h
#define Metaphor_opt_levmarq_h

namespace im
{
    // Levenberg-Mrquardt solver
    
    // Derive a class from this solver and implement the function to calculate the error residual vector
    // Call init with a pre-allocated vector to hold the value of x
    // Remember to initialize it to your starting point guess
    // Call step() in a while loop until it returns true, i.e while(!step());
    // You can print the progress by adding print statements in the while loop that examine error(), delta_error(), delta_x()
    // Get the result by calling state()

    struct LevenbergMarquardtParams
    {
        LevenbergMarquardtParams()
        {
            use_marquardt_mode = false; // enable marquardt conditioning
            fail_on_lambda_max = true; // quit if lambda hits max, else cap lambda and continue
            lambda_start = 1e-3; // initial value for lambda
            lambda_max = 1e8; // largest value of lambda before giving up
            lambda_ratio = 10.0; // lambda multiplies or divides by this ratio
            error_ratio_min = 1e-5; // smallest error vector change ratio for termination
            deriv_delta_ratio = 1e-6; // numerical deriv delta as a ratio of parameter size
            deriv_delta_min = 1e-6; // min value of deriv delta
            update_fraction = 1.0; // amount of update to add each time (<1 means don't step the whole way)
        }
        
        bool use_marquardt_mode;
        bool fail_on_lambda_max;
        double lambda_start;
        double lambda_max;
        double lambda_ratio;
        double error_ratio_min;
        double deriv_delta_ratio;
        double deriv_delta_min;
        double update_fraction;
    };
    
    // Minimizes the least squares error ||f(X)||
    // Here f(X) is a non-linear vector "residual" function of the parameter vector X.
    // You must supply the evaluator for the vector error residual f(X).
    // You can optionally supply functions for numerically or directly computing the Jacobian of f(X).
    
    template <typename TT>
    class LevenbergMarquardt
    {
    public:
        LevenbergMarquardt(LevenbergMarquardtParams const &p) : m_params(p) {}
        LevenbergMarquardt() {}
        virtual ~LevenbergMarquardt() {}
        
        // If you dont want to use the default parameters, set them before calling init()
        void set_parameters(LevenbergMarquardtParams const &p) { m_params = p; }
        
        // Initalize.
        // You can wrap external memory or pass in a freshly allocated vector
        void init(Vec<TT> vstate);
        
        // Take one minimization step. Call this in a loop until it returns true
        // Then get the result by calling state()
        bool step();
        
        // Returns true if the solver stopped before convergence.
        bool early_exit(){ return m_early_exit; }
        
        // Current sum squared error
        TT error() { return m_error; }
        
        // Get the change in error during the last step
        TT delta_error() { return m_delta_error; }
        
        // Get the distance moved during the last step
        TT delta_x() { return m_delta_x; }
        
        // Diagnostic
        TT lambda() { return m_lambda; }
        
        // Get/set the current state vector x
        Vec<TT> state() { return m_vstate; }
        
        // Dimenstionality of the state vector x
        int dims() { return m_vstate.rows(); }
        
        // Functions which are to be over-ridden by derived class.
        // The only essential function is eval_residual.
        
        // Called after init is called
        virtual void eval_init() {}
        
        // Called at the start of each step.
        virtual void eval_start_step() {}
        
        // Compute the residual vector f(X), given vector X.
        virtual void eval_residual(Vec<TT> &vresidual, Vec<TT> const &vx) = 0;
        
        // Optional function for numerically computing the Jacobian
        // This function must evaluate f(X+delta) - f(X) in order to compute one column of the Jacobian
        // f(X) is provided in the vfx parameter.
        // Typically you would add delta to vx(jcol) and then compute vdiff = f(vx) - vfx.
        // You can change vx(jcol) in this way but do not change any other elements of vx.
        virtual void eval_residual_diff(Vec<TT> &vdiff, Vec<TT> &vx, int jcol, TT delta, Vec<TT> const &vfx);
        
        // Optional function for directly computing the Jacobian
        // The default implementation computes it numerically by calling eval_residual_diff for each column.
        // The Jacobian matrix is defined as J_ij = df_i / dx_j
        virtual void eval_jacobian(Mtx<TT> &mjacobian, Vec<TT> &vx);
        
        // Called at the end of each step. Return true if you need to early exit.
        virtual bool eval_end_step() { return false; }
        
    protected:
        LevenbergMarquardtParams m_params;
        MatrixDecompLDLT<TT> m_ldlt;
        Vec<TT> m_vstate;
        Vec<TT> m_vresidual;
        Vec<TT> m_vnewresidual;
        Mtx<TT> m_mJ;
        Mtx<TT> m_mJtJ;
        Vec<TT> m_vJtJdiag;
        Vec<TT> m_vJte;
        bool m_early_exit;
        TT m_error;
        TT m_delta_error;
        TT m_delta_x;
        TT m_lambda;
    };
    
}


#endif
