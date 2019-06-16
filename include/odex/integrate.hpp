
#ifndef ODEX_INTEGRATE_HPP
#define ODEX_INTEGRATE_HPP

#include "odex/make_extrapolation_stepper.hpp"
#include <cstddef>
#include <utility>


namespace odex {

/// Integrate the differential system using sensible defaults for the 
/// extrapolation scheme.
/// \param system Time derivative operator.
/// \param state Initial state of the system.
/// \param t Initial time to evaluate the system.
/// \param dt Time step size.
/// \param n Number of time steps.
/// \param observer Observer to record output at each time step.
/// \param order Order of accuracy of the extrapolation scheme
/// \param num_cores Maximum number of cores the scheme may run on
/// \param parallel Flag to distribute work across cores
template <class Weight=double, class System, class State, class Time, class NumSteps, class Observer>
State integrate(System&& system, State const& state, Time t, Time dt, NumSteps n, Observer&& observer, 
                std::size_t order=8, std::size_t num_cores=3, bool parallel=true)
{
    auto exstepper = make_extrapolation_stepper<Weight>(std::forward<System>(system), state, order, num_cores, parallel);
    
    // copy the initial state
    State y(state);

    // run the stepper in place
    exstepper.step(y, t, dt, n, std::forward<Observer>(observer));

    // return the final output
    return y;
}

} // namespace odex

#endif // ODEX_INTEGRATE_HPP

