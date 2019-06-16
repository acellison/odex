
#ifndef ODEX_GBS_HPP
#define ODEX_GBS_HPP

#include <type_traits>
#include <cstddef>
#include <cassert>
#include <array>

namespace odex {
namespace steppers {

/// Gragg-Bulirsch-Stoer time stepper.  The asympotic error expansion contains
/// even-order terms only so that each extrapolation yields two orders of
/// accuracy.  In addition, their extrapolates have good imaginary axis 
/// coverage and are therefore useful in Method Of Lines algorithms for solving
/// hyperbolic PDE.  
template <class StateType>
class gbs
{
public:
    using state_type = StateType;
    using scratch_type = std::array<state_type, 3>;

    template <class System, class Time, class Subintervals, class SystemResult>
    static void step(System&& system, state_type const& y0, state_type& y, Time t, Time dt, Subintervals n, SystemResult&& fval0, scratch_type& scratch)
    {
        auto const h = dt/static_cast<float>(n);
        auto tn = static_cast<decltype(h)>(t);

        // Initial Forward Euler Step
        scratch[0] = y0 + h*std::forward<SystemResult>(fval0);
        tn += h;

        // First Leap Frog Step to avoid the initial data copy
        scratch[1] = y0 + 2*h*std::forward<System>(system)(tn, scratch[0]);

        // Leap Frog Iteration
        constexpr std::array<std::array<std::size_t,3>,3> inds{{ {{0, 1, 2}}, {{1, 2, 0}}, {{2, 0, 1}} }};
        std::size_t cur = 2;
        for (std::size_t ii = 1; ii < n; ++ii)
        {
            cur = cur<2 ? cur+1 : 0;
            tn += h;
            auto const ind0 = inds[cur][0];
            auto const ind1 = inds[cur][1];
            auto const ind2 = inds[cur][2];
            scratch[ind2] = scratch[ind0] + 2*h*std::forward<System>(system)(tn, scratch[ind1]);
        }

        // Smoothing Step
        auto const ind0 = inds[cur][0];
        auto const ind1 = inds[cur][1];
        auto const ind2 = inds[cur][2];
        y = .25f*(scratch[ind0]+2.0f*scratch[ind1]+scratch[ind2]);
    }
};

} // namespace steppers
} // namespace odex

#endif // ODEX_GBS_HPP

