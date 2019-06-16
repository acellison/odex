
#ifndef ODEX_EXTRAPOLATION_STEPPER_HPP
#define ODEX_EXTRAPOLATION_STEPPER_HPP

#include "odex/threading/pool.hpp"
#include "odex/detail/partition.hpp"
#include "odex/observers/null_observer.hpp"
#include <iterator>
#include <cstddef>
#include <vector>
#include <memory>

namespace odex {

/// Extrapolation stepper object.  Renders individual time stepping routines
/// at varying time step sizes, then combines the results to achieve higher
/// order accuracy by canceling terms in the asymptotic error expansions.
/// This is an underdetermined extrapolation stepper- rather than using
/// precisely enough threads to cancel that many terms in the asymptotic
/// error formula, more than necessary are utilized.  This trades overall
/// work for optimization of the stability domain of the algorithm.  Weights
/// can be computed to maximize the stability domain of the extrapolation
/// scheme over the imaginary axis for hyperbolic PDE, or over the negative
/// real axis for parabolic PDE, when used in conjuction with Method Of Lines.
/// The benefit is that each extrapolation thread can be computed completely
/// independent of the others, so on a multicore machine overall time to
/// solution is simultaneously reduced while maximizing the time step size.
template <class System, class Stepper, class State, class Weight>
class extrapolation_stepper
{
public:
    using system_type = System;
    using stepper_type = Stepper;
    using state_type = State;
    using weight_type = Weight;
    using stepper_scratch_type = typename stepper_type::scratch_type;

    /// Construct the extrapolation stepper object.
    /// \param stepper Time stepper object that does the actual system evaluation.
    /// \param system Derivative function that takes time and state.
    /// \param num_steppers Number of individual time steppers in the extrapolation scheme.
    /// \param step_counts Number of step counts for each stepper.
    /// \param weights Extrapolation weights for the output of each stepper.
    /// \param order Order of accuracy of the extrapolation scheme.
    /// \param isbn Normalized Imaginary Stability Boundary of the scheme.
    /// \param parallel Flag to distribute work across cores.
    template <class StepperType, class SystemType, class StepCountIterator, class WeightIterator>
    extrapolation_stepper(StepperType&& stepper, SystemType&& system, std::size_t num_steppers,
                          StepCountIterator step_counts, WeightIterator weights,
                          std::size_t order, float isbn, bool parallel)
    : m_order(order)
    , m_isbn(isbn)
    , m_stepper(std::forward<StepperType>(stepper))
    , m_systems()
    , m_scratch()
    , m_weights(weights, weights+static_cast<std::ptrdiff_t>(num_steppers))
    , m_step_counts(step_counts, step_counts+static_cast<std::ptrdiff_t>(num_steppers))
    , m_outputs(num_steppers)
    , m_input(nullptr)
    , m_t(0)
    , m_dt(0)
    , m_pool(nullptr)
    {
        if (parallel)
        {
            _initialize_pool(std::forward<SystemType>(system));
        }
        else
        {
            m_systems.emplace_back(std::forward<SystemType>(system));
            m_scratch.emplace_back(stepper_scratch_type());
        }
    }

    /// Order of accuracy of the time stepping scheme
    std::size_t order() const
    {
        return m_order;
    }

    /// Normalized Imaginary Stability Boundary of the scheme
    float isbn() const
    {
        return m_isbn;
    }

    /// Step the system n time steps without observation.
    /// \param y Input/output state.
    /// \param t Initial time for system evaluation.
    /// \param dt Time step size.
    /// \param n Number of time steps.
    template <class Time, class NumSteps>
    void step(state_type& y, Time t, Time dt, NumSteps n)
    {
        step(y, t, dt, n, observers::null_observer{});
    }

    /// Step the system a n time steps, observing each output.
    /// \param y Input/output state.
    /// \param t Initial time for system evaluation.
    /// \param dt Time step size.
    /// \param n Number of time steps.
    /// \param observer Callable observer object to record each time step.
    template <class Time, class NumSteps, class Observer>
    void step(state_type& y, Time t, Time dt, NumSteps n, Observer&& observer)
    {
        auto const nsteppers = m_step_counts.size();
        auto const& weights = m_weights;
        auto const& outputs = m_outputs;

        for (std::size_t ii = 0; ii < n; ++ii)
        {
            // Run the individual steppers.  This should be done depending on
            // execution context: parallel or single threaded
            _evaluate(y, t, dt);

            // Extrapolate the results from the individual steppers to get the
            // high-order-accurate result with desired stability domain.
            y = weights[0]*outputs[0];
            for (std::size_t jj = 1; jj < nsteppers; ++jj)
            {
                y += weights[jj]*outputs[jj];
            }

            // Send the result to the observer
            std::forward<Observer>(observer)(t, y);

            // Step time forward
            t = t+dt;
        }
    }

private:
    /// Dispatch the actual time stepper evaluation code.
    template <class Time>
    void _evaluate(state_type const& y, Time t, Time dt)
    {
        m_input = &y;
        m_t = t;
        m_dt = dt;
        if (m_pool)
        {
            _evaluate_parallel();
        }
        else
        {
            _evaluate_serial();
        }
    }

    /// Run the time steppers all on a single core.
    void _evaluate_serial()
    {
        // get references to data members
        auto& system = m_systems[0];
        auto const& input = *m_input;
        auto& outputs = m_outputs;
        auto const t = m_t;
        auto const dt = m_dt;
        auto const& step_counts = m_step_counts;
        auto& scratch = m_scratch[0];

        // evaluate the system to share with all steppers
        auto fval0 = system(t, input);

        // run the individual time steppers
        for (std::size_t jj = 0, nsteppers = m_step_counts.size(); jj < nsteppers; ++jj)
        {
            m_stepper.step(system, input, outputs[jj], t, dt, step_counts[jj], fval0, scratch);
        }
    }
    
    /// Run the time steppers in parallel across cores.
    void _evaluate_parallel()
    {
        m_pool->process();
    }

    /// Initialize the thread pool, dividing up the work as evenly as possible
    /// among the cores.
    template <class SystemType>
    void _initialize_pool(SystemType&& system)
    {
        // compute the core partitioning
        m_partitions = detail::partition(m_step_counts.begin(), m_step_counts.size());
        auto num_cores = m_partitions.size();

        // grab the indices corresponding to the steppers on each partition
        m_partition_indices = m_partitions;
        for (std::size_t ii = 0; ii < m_partition_indices.size(); ++ii)
        {
            auto const& partition = m_partitions[ii];
            for (std::size_t jj = 0; jj < partition.size(); ++jj)
            {
                auto iter = std::find(m_step_counts.begin(), m_step_counts.end(), partition[jj]);
                auto index = std::distance(m_step_counts.begin(), iter);
                m_partition_indices[ii][jj] = static_cast<std::size_t>(index);
            }
        }

        // target work function
        auto target = [this](std::size_t index)
        {
            // get the partition-local indices
            auto const& inds = m_partition_indices[index];

            // get local references to data members
            auto& current_system = m_systems[index];
            auto const& input = *m_input;
            auto& outputs = m_outputs;
            auto const t = m_t;
            auto const dt = m_dt;
            auto const& step_counts = m_step_counts;
            auto& scratch = m_scratch[index];

            // evaluate the system to share with all steppers on this core
            auto fval0 = current_system(t, input);

            // run each of the steppers on this core
            for (std::size_t jj = 0; jj < inds.size(); ++jj)
            {
                auto ind = inds[jj];
                m_stepper.step(current_system, input, outputs[ind], t, dt, step_counts[ind], fval0, scratch);
            }
        };

        // instantiate the thread pool
        m_pool.reset(new threading::pool(m_partitions.size()));

        // construct the workers with the target functions
        for (std::size_t ii = 0; ii < num_cores; ++ii)
        {
            m_pool->emplace(ii, target, ii);
        }

        // copy the system for each thread
        m_systems = std::vector<system_type>(m_partitions.size(), std::forward<SystemType>(system));
        m_scratch = std::vector<stepper_scratch_type>(m_partitions.size());
    }

private:
    /// order of the time stepping scheme
    std::size_t const m_order;

    /// Normalized Imaginary Stability Boundary of the scheme
    float m_isbn;

    /// time stepping algorithm.  evaluating its step() method must not change
    /// any of its internal state since this is run concurrently
    stepper_type const m_stepper;

    /// vector of systems to time step, one copy per core, so that evaluation
    /// can be performed concurrently without worrying about clobbering internal
    /// state.  if system evaluation is reentrant, consider wrapping this in
    /// a std::reference_wrapper to share this read-only memory across cores
    std::vector<system_type> m_systems;

    /// each core gets a copy of the scratch required by the time stepper
    std::vector<stepper_scratch_type> m_scratch;

    /// extrapolation weights
    std::vector<weight_type> const m_weights;

    /// extrapolation step count sequence
    std::vector<std::size_t> const m_step_counts;

    /// pre-extrapolated outputs for each time stepper
    std::vector<state_type> m_outputs;

    /// pointer to the current input
    state_type const* m_input;

    /// current time
    weight_type m_t;

    /// time step size
    weight_type m_dt;

    /// core partitions containing the step count sequence for each core
    std::vector<std::vector<std::size_t>> m_partitions;

    /// indices into the weights and step counts of each time stepper.
    /// ordering matches that of the m_partitions vector to enable easy
    /// lookup into the linear arrays
    std::vector<std::vector<std::size_t>> m_partition_indices;

    /// thread pool that dispatches the workers
    std::unique_ptr<threading::pool> m_pool;
};

} // namespace odex

#endif // ODEX_EXTRAPOLATION_STEPPER_HPP
