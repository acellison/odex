
#ifndef ODEX_OBSERVERS_NULL_OBSERVER_HPP
#define ODEX_OBSERVERS_NULL_OBSERVER_HPP

namespace odex {
namespace observers {

/// No-op observer
struct null_observer
{
    template <class Time, class System>
    constexpr void operator()(Time&&, System&&) { }
};

} // namespace observers
} // namespace odex

#endif // ODEX_OBSERVERS_NULL_OBSERVER_HPP
