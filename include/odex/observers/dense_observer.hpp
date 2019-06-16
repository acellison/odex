
#ifndef ODEX_OBSERVERS_DENSE_OBSERVER_HPP
#define ODEX_OBSERVERS_DENSE_OBSERVER_HPP

#include <vector>

namespace odex {
namespace observers {

/// A dense observer records the time and state after each sample
/// is computed by the extrapolation_stepper.
template <class Time, class State>
class dense_observer
{
public:
    dense_observer()
    {    }

    explicit dense_observer(std::size_t size)
    {
        m_time.reserve(size);
        m_state.reserve(size);
    }

    template <class T, class S>
    void operator()(T&& t, S&& state)
    {
        m_time.push_back(std::forward<T>(t));
        m_state.push_back(std::forward<S>(state));
    }

    std::vector<T> const& time() const { return m_time; }
    std::vector<T>      & time()       { return m_time; }
    std::vector<T> const& state() const { return m_state; }
    std::vector<T>      & state()       { return m_state; }

private:
    std::vector<Time> m_time;
    std::vector<State> m_state;
};

} // namespace observers
} // namespace odex

#endif // ODEX_OBSERVERS_DENSE_OBSERVER_HPP

