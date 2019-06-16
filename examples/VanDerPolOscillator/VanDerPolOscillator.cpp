#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Core>

using ValueType = double;
using StateType = Eigen::Array<ValueType,2,1>;

class VanDerPolOscillator
{
public:
    explicit VanDerPolOscillator(ValueType mu)
    : m_mu(mu)
    {

    }

    StateType operator()(ValueType, StateType const& state) const
    {
        ValueType x = state[0];
        ValueType y = state[1];
        StateType yp;
        yp[0] = y;
        yp[1] = m_mu*(1-x*x)*y-x;
        return yp;
    }

private:
    ValueType m_mu;
};

int main()
{
    // Set up the ODE
    ValueType mu = 10.65;
    VanDerPolOscillator system(mu);

    // Stepper parameters
    std::size_t nsteps = 16*2048;
    double t0 = 0;
    double t1 = 100;
    double dt = (t1-t0)/nsteps;

    // Initial state
    StateType y0; y0[0] = 1; y0[1] = 0;

    // Observer records each state
    std::vector<ValueType> tn(nsteps), xn(nsteps), yn(nsteps);
    std::size_t index = 0;
    auto observer = [&index, &tn, &xn, &yn](auto t, auto const& y)
    {
        tn[index] = t;
        xn[index] = y[0];
        yn[index] = y[1];
        ++index;
    };

    // Run the odex numerical integration
    odex::integrate(system, y0, t0, dt, nsteps, observer);

    // Plot
    namespace plt = matplotlibcpp;
    plt::backend("Qt5Agg");

    plt::figure();
    plt::plot(tn, xn);
    plt::plot(tn, yn);
    plt::xlabel("Time");
    plt::ylabel("Value");
    plt::title("Van Der Pol Oscillator Time Series");
    plt::legend();
    plt::grid(true);

    plt::figure();
    plt::plot(xn, yn);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::title("Van Der Pol Oscillator Phase Portrait");
    plt::grid(true);
    plt::show();
}
