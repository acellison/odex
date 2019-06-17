#include "odex/integrate.hpp"
#include "matplotlibcpp.h"
#include <Eigen/Dense>
#include <exception>
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

using ValueType = double;
using StateType = Eigen::Array<ValueType,1000,1>;

class KdV
{
public:
    explicit KdV(ValueType c, ValueType k, ValueType xmin, ValueType xmax, ValueType xi0)
    : m_c(c), m_k(k), m_xmin(xmin), m_xmax(xmax), m_xi0(xi0)
    {

    }

    auto operator()(ValueType t, StateType const& u)
    {
        deriv1_4(t, u, m_ux);
        deriv3_4(t, u, m_uxxx);
        return -(6*u*m_ux + m_uxxx);
    }

    /// Fourth order accurate centered finite difference
    /// approximation to the first derivatve
    void deriv1_4(ValueType t, StateType const& u, StateType& ux)
    {
        auto n = u.size();
        auto k = m_k;
        auto unp2 = u.segment(4,n-4);
        auto unp1 = u.segment(3,n-4);
        auto unm1 = u.segment(1,n-4);
        auto unm2 = u.segment(0,n-4);
        ux.segment(2,n-4) = (unm2-8*unm1+8*unp1-unp2)/(12*k);

        // boundary values are the true solution
        ValueType un1 = soliton_value(m_c, m_xmin-1*k-m_c*t-m_xi0);
        ValueType un2 = soliton_value(m_c, m_xmin-2*k-m_c*t-m_xi0);
        ValueType up1 = soliton_value(m_c, m_xmax+1*k-m_c*t-m_xi0);
        ValueType up2 = soliton_value(m_c, m_xmax+2*k-m_c*t-m_xi0);
        ux[1]   = (un1   -8*u[0]  +8*u[2]  -u[3])/(12*k);
        ux[0]   = (un2   -8*un1   +8*u[1]  -u[2])/(12*k);
        ux[n-1] = (u[n-3]-8*u[n-2]+8*up1   -up2 )/(12*k);
        ux[n-2] = (u[n-4]-8*u[n-3]+8*u[n-1]-up1 )/(12*k);
    }

    /// Fourth order accurate centered finite difference
    /// approximation to the third derivatve
    void deriv3_4(ValueType t, StateType const& u, StateType& uxxx)
    {
        auto n = u.size();
        auto k = m_k;
        auto k3 = m_k*m_k*m_k;
        auto unp3 = u.segment(6,n-6);
        auto unp2 = u.segment(5,n-6);
        auto unp1 = u.segment(4,n-6);
        auto unm1 = u.segment(2,n-6);
        auto unm2 = u.segment(1,n-6);
        auto unm3 = u.segment(0,n-6);
        uxxx.segment(3,n-6) = (unm3-8*unm2+13*unm1-13*unp1+8*unp2-unp3)/(8*k3);

        // boundary values are the true solution
        ValueType un1 = soliton_value(m_c, m_xmin-1*k-m_c*t-m_xi0);
        ValueType un2 = soliton_value(m_c, m_xmin-2*k-m_c*t-m_xi0);
        ValueType un3 = soliton_value(m_c, m_xmin-3*k-m_c*t-m_xi0);
        ValueType up1 = soliton_value(m_c, m_xmax+1*k-m_c*t-m_xi0);
        ValueType up2 = soliton_value(m_c, m_xmax+2*k-m_c*t-m_xi0);
        ValueType up3 = soliton_value(m_c, m_xmax+3*k-m_c*t-m_xi0);
        uxxx[2]   = (un1   -8*u[0]  +13*u[1]  -13*u[3]  +8*u[4]  -u[5])/(8*k3);
        uxxx[1]   = (un2   -8*un1   +13*u[0]  -13*u[2]  +8*u[3]  -u[4])/(8*k3);
        uxxx[0]   = (un3   -8*un2   +13*un1   -13*u[1]  +8*u[2]  -u[3])/(8*k3);
        uxxx[n-1] = (u[n-4]-8*u[n-3]+13*u[n-2]-13*up1   +8*up2   -up3 )/(8*k3);
        uxxx[n-2] = (u[n-5]-8*u[n-4]+13*u[n-3]-13*u[n-1]+8*up1   -up2 )/(8*k3);
        uxxx[n-3] = (u[n-6]-8*u[n-5]+13*u[n-4]-13*u[n-2]+8*u[n-1]-up1 )/(8*k3);
    }

    /// Second order accurate centered finite difference
    /// approximation to the first derivatve
    void deriv1_2(StateType const& u, StateType& ux)
    {
        auto n = u.size();
        ux.segment(1,n-2) = (u.tail(n-2)-u.head(n-2))/(2*m_k);
        ux[0]   = (u[1]       )/(2*m_k);
        ux[n-1] = (    -u[n-2])/(2*m_k);
    }

    /// Second order accurate centered finite difference
    /// approximation to the third derivatve
    void deriv3_2(StateType const& u, StateType& uxxx)
    {
        auto n = u.size();
        auto k3 = m_k*m_k*m_k;
        auto unp2 = u.segment(4,n-4);
        auto unp1 = u.segment(3,n-4);
        auto unm1 = u.segment(1,n-4);
        auto unm2 = u.segment(0,n-4);
        uxxx.segment(2,n-4) = (unm2-2*unm1+2*unp1-unp2)/(-2*k3);
        uxxx[1]   = (      -2*u[0]  +2*u[2]  -u[3])/(-2*k3);
        uxxx[0]   = (                2*u[1]  -u[2])/(-2*k3);
        uxxx[n-1] = (u[n-3]-2*u[n-2]              )/(-2*k3);
        uxxx[n-2] = (u[n-4]-2*u[n-3]+2*u[n-1]     )/(-2*k3);
    }

    static ValueType soliton_value(ValueType c, ValueType xi)
    {
        auto s = 1/std::cosh(std::sqrt(c)/2*xi);
        return c*(s*s)/2;
    }

private:
    ValueType const m_c;
    ValueType const m_k;
    ValueType const m_xmin;
    ValueType const m_xmax;
    ValueType const m_xi0;
    StateType m_ux;
    StateType m_uxxx;
};

int main()
{
    // odex integrator parameters
    std::size_t order = 8;
    std::size_t cores = 3;

    // Set up the PDE
    ValueType c = 1;
    ValueType xmin = -100;
    ValueType xmax =  100;
    ValueType xi0 = xmin/2;
    ValueType xifinal = xmax/2;
    auto npoints = StateType::SizeAtCompileTime;
    auto k = (xmax-xmin)/(npoints-1);

    // Stepper parameters
    double t0 = 0;
    double dt = 1e-2;
    std::size_t nsteps = std::size_t(std::round((xifinal-xi0)/(c*dt)));

    // Construct the soliton with initial displacement a
    auto soliton = [c,k,xmin](ValueType a)
    {
        StateType u;
        auto n = u.size();
        for (decltype(n) ii = 0; ii < n; ++ii)
        {
            auto xi = xmin+k*ii-a;
            u[ii] = KdV::soliton_value(c, xi);
        }
        return u;
    };

    // Initial state
    StateType u0 = soliton(xi0);

    // Construct the system
    KdV system(c, k, xmin, xmax, xi0);

    // Observer records a subset of the actual outputs
    std::size_t decfactor = 1e3;
    std::vector<std::vector<double>> un(nsteps/decfactor+1, std::vector<double>(std::size_t(npoints), 0));;
    std::copy_n(u0.data(), u0.size(), un[0].begin());
    std::size_t index = 0;
    auto observer = [&index, &un, decfactor](auto, auto const& u)
    {
        ++index;
        if (index % decfactor == 0)
        {
            if (std::any_of(u.data(), u.data()+u.size(),
                            [](ValueType v){ return std::isnan(v) || std::isinf(v); }))
            {
                throw std::runtime_error("Solution is unstable!");
            }
            std::copy_n(u.data(), u.size(), un[index/decfactor].begin());
        }
    };

    // Run the odex numerical integration
    StateType ufinal = odex::integrate<ValueType>(system, u0, t0, dt, nsteps, observer, order, cores);

    // Compute the error relative to the true solution
    auto norm = [](StateType a)
    {
        return a.matrix().lpNorm<Eigen::Infinity>();
    };
    StateType ufinal_true = soliton(xifinal);
    StateType error = ufinal_true-ufinal;
    auto linf_error = norm(error)/0.5;  // normalize to max of soliton
    std::cout << "Relative Soliton L-Inf Error: " << linf_error << std::endl;

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
    plt::title("KdV Soliton Solution");
    plt::grid(true);

    plt::figure();
    std::vector<double> verror(error.data(), error.data()+npoints);
    plt::plot(verror);
    plt::title("Error in Soliton Solution");
    plt::xlabel("x");
    plt::ylabel("error");
    plt::grid(true);

    plt::show();
}
