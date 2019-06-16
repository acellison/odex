#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Core>

using ValueType = double;
using StateType = Eigen::Matrix<ValueType,2,256,Eigen::RowMajorBit>;

class WaveEquation
{
public:
    explicit WaveEquation(ValueType c, ValueType k)
    : m_c2(c*c), m_k(k)
    {

    }

    auto const& operator()(ValueType, StateType const& u)
    {
        m_state.row(0) = u.row(1);
        gradient2(u.row(0), m_state.row(1));
        m_state.row(1) *= m_c2;
        return m_state;
    }

    template <class Lhs, class Rhs>
    void gradient2(Lhs const& u, Rhs&& uxx)
    {
        auto n = m_state.cols();
        auto scale = 1/(m_k*m_k);
        uxx.segment(1,n-2) = (u.tail(n-2)-2*u.segment(1,n-2)+u.head(n-2))*scale;

        // Zero-displacement boundary conditions
        uxx[0]   = (u[1]-2*u[0]  +0     )*scale;
        uxx[n-1] = (0   -2*u[n-1]+u[n-2])*scale;
    }

private:
    ValueType m_c2;
    ValueType m_k;
    StateType m_state;
};

int main()
{
    // Set up the PDE
    ValueType c = 1.0;
    ValueType k = 1.0;
    WaveEquation system(c,k);

    // Stepper parameters
    std::size_t nsteps = 16384;
    double t0 = 0;
    double t1 = 512;
    double dt = (t1-t0)/nsteps;

    // Initial state
    StateType u0 = StateType::Zero();
    auto n = u0.cols();
    for (decltype(n) ii = 0; ii < n; ++ii)
    {
        ValueType x = ValueType(ii)/n-.5;
        u0(0,ii) = std::exp(-1200*(x*x));
    }

    // Observer records a subset of the actual outputs
    std::size_t decfactor = 128;
    std::vector<std::vector<ValueType>> un(nsteps/decfactor, std::vector<ValueType>(std::size_t(n), 0));;
    std::size_t index = 0;
    auto observer = [&index, &un, decfactor](auto, auto const& u)
    {
        auto ncols = u.cols();
        if (index % decfactor == 0)
        {
            std::copy_n(u.row(0).data(), ncols, un[index/decfactor].begin());
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
    plt::title("Wave Equation Solution");
    plt::grid(true);
    plt::show();
}
