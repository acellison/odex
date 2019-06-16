#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Core>
#include <map>

using ValueType = double;
using StateType = Eigen::Array<ValueType,3,1>;

class LorenzAttractor
{
public:
    LorenzAttractor(ValueType sigma, ValueType rho, ValueType beta)
    : m_sigma(sigma), m_rho(rho), m_beta(beta)
    {

    }

    StateType operator()(ValueType, StateType const& state) const
    {
        ValueType x = state[0];
        ValueType y = state[1];
        ValueType z = state[2];
        StateType yp;
        yp[0] = m_sigma*(y-x);
        yp[1] = x*(m_rho-z)-y;
        yp[2] = x*y-m_beta*z;
        return yp;
    }

private:
    ValueType m_sigma;
    ValueType m_rho;
    ValueType m_beta;
};

int main()
{
    // Set up the ODE
    ValueType sigma = 10.0;
    ValueType rho   = 28.0;
    ValueType beta  = 8.0/3;
    LorenzAttractor system(sigma, rho, beta);

    // Stepper parameters
    std::size_t nsteps = 10000;
    double t0 = 0;
    double t1 = 100;
    double dt = (t1-t0)/nsteps;

    // Initial state
    StateType y0; y0[0] = 1; y0[1] = 0; y0[2] = 0;

    // Observer records the state at each time step
    std::vector<ValueType> xn(nsteps), yn(nsteps), zn(nsteps);
    std::size_t index = 0;
    auto observer = [&index, &xn, &yn, &zn](auto, auto const& y)
    {
        xn[index] = y[0];
        yn[index] = y[1];
        zn[index] = y[2];
        ++index;
    };

    // Run the odex numerical integration
    odex::integrate(std::ref(system), y0, t0, dt, nsteps, observer);

    // Plot
    namespace plt = matplotlibcpp;
    plt::backend("Qt5Agg");

    std::map<std::string, std::string> kwargs = {
        { "linewidth", "0.3" }
    };
    plt::plot(xn, yn, zn, kwargs);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::zlabel("z");
    plt::title("Lorenz Attractor");
    plt::grid(true);
    plt::show();
}
