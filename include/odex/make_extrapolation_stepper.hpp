
#ifndef ODEX_MAKE_EXTRAPOLATION_STEPPER_HPP
#define ODEX_MAKE_EXTRAPOLATION_STEPPER_HPP

#include "odex/extrapolation_stepper.hpp"
#include "odex/steppers/gbs.hpp"
#include "odex/detail/make_extrap_config.hpp"
#include <type_traits>
#include <utility>
#include <cstddef>
#include <vector>
#include <tuple>

namespace odex {


/// Construct an extrapolation_stepper with the given system and state.  The
/// num_cores parameter is used when selecting the extrapolation scheme's
/// weights.  If parallel is false the algorithm will run on a single core,
/// but will still use the weights resulting from the order/num_cores
/// combination.  Higher number of cores yields higher ISBn so larger time
/// steps can be taken when solving a wave-type PDE with method-of-lines.
/// \param system Time derivative operator.
/// \param state Initial state of the system.
/// \param order Order of accuracy of the extrapolation scheme.
/// \param num_cores Maximum number of cores the scheme may run on.
/// \param parallel Flag to distribute work across cores.
template <class Weight=double, class System, class State>
auto make_extrapolation_stepper(System&& system, State const& state, std::size_t order=8, std::size_t num_cores=3, bool parallel=true)
{
    using weight_type = Weight;
    using stepper_type = odex::steppers::gbs<State>;
    using exstepper_type = odex::extrapolation_stepper<std::decay_t<System>, stepper_type, State, weight_type>;

    // avoid unused parameter warning
    (void)state;

    // get the extrapolation configuration for the specified order and number of cores
    float isbn = 0.0f;
    std::vector<std::size_t> step_counts;
    std::vector<weight_type> weights;
    std::tie(isbn, step_counts, weights) = detail::make_extrap_config<weight_type>(order, num_cores);

    // construct the extrapolation stepper
    return exstepper_type(stepper_type(), std::forward<System>(system), step_counts.size(), step_counts.begin(), weights.begin(),
                          order, isbn, parallel);
}


} // end namespace odex

#endif // ODEX_MAKE_EXTRAPOLATION_STEPPER_HPP

