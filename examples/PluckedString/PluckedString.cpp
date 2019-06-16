#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Core>

using ValueType = double;
constexpr std::size_t npoints = 256;
using StateType = Eigen::Matrix<ValueType,2,npoints,Eigen::RowMajorBit>;

class PluckedString
{
public:
    explicit PluckedString(ValueType c, ValueType b, ValueType k)
    : m_c2(c*c), m_b(b), m_k(k)
    {

    }

    auto const& operator()(ValueType, StateType const& u)
    {
        m_state.row(0) = u.row(1);
        gradient2(u.row(0), m_state.row(1));
        m_state.row(1) *= m_c2;
        m_state.row(1) -= m_b*m_state.row(0);
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
    ValueType m_b;
    ValueType m_k;
    StateType m_state;
};

int main()
{
    using Index = typename StateType::Index;

    // Set up the PDE
    ValueType c = 1.0;
    ValueType b = .025;
    ValueType length = 1.0;
    ValueType k = length/npoints;
    PluckedString system(c,b,k);

    // Stepper parameters
    std::size_t nsteps = 16384;
    double t0 = 0;
    double t1 = 2;
    double dt = (t1-t0)/nsteps;

    // Initial state
    StateType u0 = StateType::Zero();
    for (std::size_t ii = 0; ii < npoints; ++ii)
    {
        ValueType x = k*ii;

        // Triangle wave
        ValueType mod = 1;
        for (std::size_t nn = 0; nn < 8; ++nn)
        {
            std::size_t m = 2*nn+1;
            u0(0,Index(ii)) += mod/(m*m)*std::sin(m*2*M_PI*x/length);
            mod *= -1;
        }
        u0(0,Index(ii)) *= 8/(M_PI*M_PI);

        // Window
        auto window = 0.5*(1+std::cos(M_PI*x/length));
        u0(0,Index(ii)) *= window*window;
    }

    // Observer records a subset of the actual outputs
    std::size_t decfactor = 256;
    std::vector<std::vector<ValueType>> un(nsteps/decfactor, std::vector<ValueType>(npoints, 0));;
    std::size_t index = 0;
    auto observer = [&index, &un, decfactor](auto, auto const& u)
    {
        auto n = u.cols();
        if (index % decfactor == 0)
        {
            std::copy_n(u.row(0).data(), n, un[index/decfactor].begin());
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
