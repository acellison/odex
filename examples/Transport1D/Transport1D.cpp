#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Core>

using ValueType = double;
using StateType = Eigen::Array<ValueType,256,1>;

class Transport1D
{
public:
    explicit Transport1D(ValueType c, ValueType k)
    : m_c(c), m_k(k)
    {

    }

    auto operator()(ValueType, StateType const& u)
    {
        auto n = u.size();
        m_ux.segment(1,n-2) = (u.tail(n-2)-u.head(n-2))/(2*m_k);
        m_ux[0]   = (u[1]-u[n-1])/(2*m_k);
        m_ux[n-1] = (u[0]-u[n-2])/(2*m_k);

        return -m_c*m_ux;
    }

private:
    ValueType m_c;
    ValueType m_k;
    StateType m_ux;
};

int main()
{
    // Set up the PDE
    ValueType c = 1.0;
    ValueType k = 1.0;
    Transport1D system(c,k);

    // Stepper parameters
    std::size_t nsteps = 2048;
    double t0 = 0;
    double t1 = 512;
    double dt = (t1-t0)/nsteps;

    // Initial state
    StateType u0;
    auto n = u0.size();
    for (decltype(n) ii = 0; ii < n; ++ii)
    {
        ValueType x = ValueType(ii)/n-.5;
        u0[ii] = std::exp(-60*(x*x));
    }

    // Observer records a subset of the actual outputs
    std::size_t decfactor = 16;
    std::vector<std::vector<ValueType>> un(nsteps/decfactor, std::vector<ValueType>(std::size_t(n), 0));;
    std::size_t index = 0;
    auto observer = [&index, &un, decfactor](auto, auto const& u)
    {
        if (index % decfactor == 0)
        {
            std::copy_n(u.data(), u.size(), un[index/decfactor].begin());
        }
        ++index;
    };

    // Run the odex numerical integration
    odex::integrate(system, u0, t0, dt, nsteps, observer);

    // Plot
    namespace plt = matplotlibcpp;
    plt::backend("Qt5Agg");

    std::map<std::string, std::string> kwargs = {
        { "linewidth", "0.8" },
        { "color",     "#1f77b4" }
    };
    plt::plot_waterfall(un, kwargs);
    plt::xlabel("x");
    plt::ylabel("Time");
    plt::zlabel("Amplitude");
    plt::title("1D Transport Solution");
    plt::grid(true);
    plt::show();
}
