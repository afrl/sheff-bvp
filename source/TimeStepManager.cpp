#include "fea_analysis/TimeStepManager.h"

namespace fea_analysis
{
    TimeStepManager::TimeStepManager()
    {
        current_time_ = 0.0;
        current_time_increment_ = 0.0;
        current_time_step_ = 0;
        previous_time_ = 0.0;
        previous_time_increment_ = 0.0;
        initial_guess_type_ = InitialGuessType::Linear;
    }

    void TimeStepManager::start_iteration_loop()
    {
        curr_iterations_.push(0);
    }

    int TimeStepManager::converged_iteration_loop()
    {
        int top = get_iteration_number();
        curr_iterations_.pop();
        return top;
    }

    bool TimeStepManager::is_first_iteration() const
    {
        return curr_iterations_.top() == 0;
    }

    int TimeStepManager::get_iteration_number() const
    {
        if (curr_iterations_.empty())
            return -1;
        else
            return curr_iterations_.top();
    }

    ///////////////////////////////////
    // StaticAnalysisTimeStepManager //
    ///////////////////////////////////

    StaticAnalysisTimeStepManager::StaticAnalysisTimeStepManager() : TimeStepManager()
    {
        current_time_step_ = 0;
        // A really big number so that an rate dependency becomes numerically zero
        current_time_increment_ = 1e10;
        max_iterations_ = 400;
    }

    void StaticAnalysisTimeStepManager::start_time_step_loop()
    {
        //update_time_scaling_curves();
        while (!curr_iterations_.empty())
            curr_iterations_.pop();
    }

    bool StaticAnalysisTimeStepManager::next_time_step()
    {
        return false;
    }

    bool StaticAnalysisTimeStepManager::next_iteration()
    {
        if (curr_iterations_.top() + 1 < max_iterations_) {
            curr_iterations_.top()++;
            return true;
        }
        else
            return false;
    }
}
