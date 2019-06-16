#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Core>

using ValueType = double;
using StateType = Eigen::Array<ValueType,2048,1>;

class ViscidBurgers
{
public:
    explicit ViscidBurgers(ValueType gamma, ValueType k)
    : m_gamma(gamma), m_k(k)
    {

    }

    auto operator()(ValueType, StateType const& u)
    {
        gradient1(u, m_ux);
        gradient2(u, m_ux, m_uxx);
        return m_gamma*m_uxx - u*m_ux;
    }

    void gradient1(StateType const& u, StateType& ux)
    {
        auto n = u.size();
        ux.segment(1,n-2) = (u.tail(n-2)-u.head(n-2))/(2*m_k);
        ux[0]   = (u[1]  -u[0]  )/(m_k);
        ux[n-1] = (u[n-1]-u[n-2])/(m_k);
    }

    void gradient2(StateType const& u, StateType const& ux, StateType& uxx)
    {
        auto n = u.size();
        uxx.segment(1,n-2) = (u.tail(n-2)-2*u.segment(1,n-2)+u.head(n-2))/(m_k*m_k);
        uxx[0]   = (ux[1]  -ux[0]  )/m_k;
        uxx[n-1] = (ux[n-1]-ux[n-2])/m_k;
    }

private:
    ValueType m_gamma;
    ValueType m_k;
    StateType m_ux;
    StateType m_uxx;
};

int main()
{
    // Set up the PDE
    ValueType gamma = 4.0;
    ValueType k = 1e-1;
    ViscidBurgers system(gamma,k);

    // Stepper parameters
    std::size_t nsteps = 1e5;
    double t0 = 0;
    double dt = 2.5e-3;

    // Initial state
    StateType u0;
    auto n = u0.size();
    for (decltype(n) ii = 0; ii < n; ++ii)
    {
        ValueType x = ValueType(ii)-n/6;
        u0[ii] = 1-std::tanh(k*x/(2*gamma));
    }

    // Observer records a subset of the actual outputs
    std::size_t decfactor = 1000;
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
    plt::title("Viscid Burgers Traveling Front Solution");
    plt::grid(true);
    plt::show();
}
