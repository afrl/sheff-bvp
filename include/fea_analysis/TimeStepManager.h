#pragma once

#include <stack>

namespace fea_analysis
{

    /// A class to manage the time step loops and iteration loops during an analysis.
    /// It should keep the state (what is the current time and iteration) inform the
    /// model if it is done time stepping or if too many iterations have passed.
    /// Derived classes should be implemented for all the different time stepping
    /// schemes.  Even static analysis will use this class to manage iteration loops,
    /// though it there is no time stepping persay.
    class TimeStepManager {
    public:
        /// Enumerates the type of intial guesses for a time step
        enum class InitialGuessType : short {
            Linear, ///< Use a linear solution as the initail guess
            Extrapolation ///< Use a linear extrapolation (same solution increment as last step as initial guess
        };

        TimeStepManager();

        //////////////////
        // Manipulators //
        //////////////////

        /// Initialize any values needed to start a new time loop.  Do not assume
        /// that this is a new object.  Reset state for a new time loop.
        virtual void start_time_step_loop() = 0;
        /// Calculates the time for the next time step.  This should be called at
        /// the end of a time loop. Returns true if there is another step left.
        /// Suggested structure: \code{cpp}
        /// ts_anager.start_time_step_loop();
        /// do {
        ///		// Apply loads, solve, write output, etc.
        /// } while (ts_anager.next_time_step())
        /// \endcode
        virtual bool next_time_step() = 0;
        /// Should be called before an iteration loop is entered.  It resets
        /// the current iteration and convergence conditions.
        /// It is fine to have nested iteration loops.  Just make sure to call
        /// TimeStepManager::converged_iteration_loop when it converges
        virtual void start_iteration_loop();
        /// Starts the next iteration.  Returns true if there is another iteration
        /// allowed according to the time stepping scheme.
        virtual bool next_iteration() = 0;
        /// Indicate that the current iteration loop is converged
        virtual int converged_iteration_loop();
        /// Attempts to cutback on the load but decreasing the time increment
        /// by a factor.
        /// \param factor - the factor for the load cutback, should be between 0 and 1
        /// Returns whether the cutback was allowed.
        virtual bool try_load_cutback(double factor) { return false; };

        ///////////////
        // Accessors //
        ///////////////

        /// Returns a bool indicating if this is the first time step of the loop.
        virtual bool is_first_time_step() const { return current_time_step_ == 0; }
        virtual bool is_last_time_step() const { return true; }
        /// Returns the current time value (or psuedo-time if it is quasi-static)
        virtual double get_current_time() const { return current_time_; }
        virtual double get_current_time_increment() const { return current_time_increment_; }
        /// Returns the current time step number
        virtual int get_time_step() const { return current_time_step_; }
        virtual double get_previous_time() { return previous_time_; }
        virtual double get_previous_time_increment() { return previous_time_increment_; }
        virtual double get_extrapolation_factor() const { return 1.0; };
        /// Returns a bool indicating if this is the first time step of the loop.
        virtual bool is_first_iteration() const;
        /// Returns the current iteration number of the current iteration loop,
        /// which would be the most nested loop.
        int get_iteration_number() const;
        inline InitialGuessType get_intial_guess_type() const { return initial_guess_type_; }

    protected:

        /// Current time of the analysis
        double current_time_;
        /// Current time increment (takes into account cutbacks)
        double current_time_increment_;
        /// The current time step that the analysis is on.
        int current_time_step_;
        /// The time at the previously converged time step
        double previous_time_;
        /// Time increment from last time step (takes into account cutbacks)
        double previous_time_increment_;
        /// A stack to hold the current iteration of all the iteration
        /// loops ongoing.  The top of the stack is the current
        /// iteration loop.
        std::stack<int> curr_iterations_;

        InitialGuessType initial_guess_type_;
    };

    /// A TimeStepManager for a StaticAnalysis.  It may sound counterintuitive since
    /// a static analysis does not have 'time' but it allows manages iterating nonlinearities
    /// and provides a consistent interface for the model to use (making the same code
    /// for quasi-static to work with static).
    /// Note: for now the max iterations is stored here, but I think there should be
    /// TimeConvergenceCriteria and IterationConvergenceCriteria classes implemented that
    /// that allow a much more flexible mix and match of convergence criteria with time
    /// stepping schemes.
    class StaticAnalysisTimeStepManager : public TimeStepManager
    {
    public:
        /// Constructor
        StaticAnalysisTimeStepManager();

        /// Restarts any iteration loops.
        virtual void start_time_step_loop();
        /// Just returns false since there is only one step for a static analysis.
        virtual bool next_time_step() override;
        /// See TimeStepManager::next_iteration
        virtual bool next_iteration() override;

    protected:
        /// The maximum number of iterations allowed for any iteration loop.
        int max_iterations_;
    };

}
