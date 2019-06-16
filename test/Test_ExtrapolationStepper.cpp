
#ifdef NDEBUG
#  undef NDEBUG
#endif // NDEBUG

#include "odex/integrate.hpp"
#include "odex/extrapolation_stepper.hpp"
#include "odex/steppers/gbs.hpp"
#include "convector.hpp"
#include "matrix.hpp"
#include <iostream>
#include <cassert>
#include <vector>
#include <array>
#include <cmath>

static void run_simple_ode(std::size_t order, std::size_t num_cores, bool parallel, bool print)
{
    using state_type = double;
    using time_type = double;

    auto system = [](auto, auto y)
    {
        return y;
    };

    std::size_t nsteps = 32;
    time_type t0  = 0;
    time_type t1  = 2;
    time_type dt  = (t1-t0)/nsteps;

    std::vector<state_type> output(nsteps+1, state_type{});
    state_type y0 = std::exp(t0);
    state_type y  = odex::integrate(system, y0, t0, dt, nsteps, odex::observers::null_observer{},
                                    order, num_cores, parallel);

    auto error = std::exp(t1)-y;
    if (print)
    {
        std::cout << "odex " << (parallel ? "parallel" : "serial  ") <<  ": simple ode error  " << error << std::endl;
    }
    assert(error < 3e-12 && "odex error too large!");
}

static void test_simple_ode()
{
    std::vector<std::array<std::size_t,2>> configs = { {8,3}, {8,6}, {8,8}, {12,4}, {12,8}, {16,5} };

    // Test each configuration
    for (std::size_t ii = 0; ii < configs.size(); ++ii)
    {
        auto order = configs[ii][0];
        auto cores = configs[ii][1];

        run_simple_ode(order, cores, false, true);
        run_simple_ode(order, cores, true,  true);
    }


    // Bang on the threading synchronization
    std::size_t order = 8;
    std::size_t cores = 8;
    bool parallel = true;
    bool print = false;
    for (std::size_t jj = 0; jj < 40; ++jj)
    {
        run_simple_ode(order, cores, parallel, print);
    }
}

static double run_convection_2d(std::size_t order, std::size_t cores, bool parallel)
{
    std::cout << "Running GBS_{" << order << "," << cores << "}: 2D Convection in " << (parallel ? "Parallel" : "Series") << "..." << std::endl;

    constexpr std::size_t npoints = 256;
    using weight_type = double;
    using value_type = double;
    using state_type = matrix<value_type, npoints, npoints>;
    using time_type = value_type;
    using system_type = convector<value_type, state_type>;

    std::size_t nsteps = 256;
    time_type t0 = 0;
    time_type t1 = 1e-3;
    time_type dt = (t1-t0)/nsteps;

    value_type k = 1;
    value_type c[] = { 0.5, 0.25 };

    system_type system(k, c[0], c[1]);

    state_type u0;
    for (std::ptrdiff_t ii = 0; ii < ptrdiff_t(npoints); ++ii)
    {
        for (std::ptrdiff_t jj = 0; jj < ptrdiff_t(npoints); ++jj)
        {
            double x = static_cast<double>(ii)/npoints-.5;
            double y = static_cast<double>(jj)/npoints-.5;
            double norm = x*x+y*y;
            u0(ii,jj) = std::exp(-60*norm);
        }
    }

    auto now = []()
    {
        return uint64_t(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    };

    // construct the extrapolation stepper
    auto exstepper = odex::make_extrapolation_stepper<weight_type>(system, u0, order, cores, parallel);

    // copy initial state
    auto u = u0;

    // run the stepper
    auto begin_time = now();
    exstepper.step(u, t0, dt, nsteps);
    auto end_time = now();
    auto duration = 1e-9*(end_time-begin_time);
    return duration;
}

static void test_convection_2d()
{
    std::vector<std::array<std::size_t,2>> configs = { {8,3}, {8,6}, {8,8}, {12,4}, {12,8}, {16,5} };

    for (std::size_t ii = 0; ii < configs.size(); ++ii)
    {
        auto order = configs[ii][0];
        auto cores = configs[ii][1];
        auto duration_serial = run_convection_2d(order,cores,false);
        auto duration_parallel = run_convection_2d(order,cores,true);
        auto speedup = duration_serial/duration_parallel;
        std::cout << "  parallel speedup: " << speedup << std::endl;
        std::cout << "  parallel efficiency: " << speedup/cores*100 << "%" << std::endl;
    }
}

int main()
{
    test_simple_ode();
    test_convection_2d();
}
