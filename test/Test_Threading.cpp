
#include "odex/threading/pool.hpp"
#include "odex/threading/worker.hpp"
#include <iostream>
#include <cassert>
#include <thread>
#include <atomic>

struct target
{
    target()
    : counter(0)
    , done(false)
    {    }

    target(target const& other)
    : counter(other.counter.load())
    , done(other.done.load())
    {    }

    void operator()()
    {
        ++counter;
        done = true;
    }

    std::atomic<int> counter;
    std::atomic<bool> done;
};

static void test_worker_notify()
{
    std::size_t const iters = 10;

    target targ;

    odex::threading::worker worker(std::ref(targ));
    for (std::size_t ii = 0; ii < iters; ++ii)
    {
        worker.notify();
        while (!targ.done)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        targ.done = false;
    }
    worker.join();

    assert(targ.counter == iters);
}

static void test_thread_pool()
{
    std::size_t const num_workers = 4;
    std::size_t const iters = 10;

    target targets[num_workers];

    odex::threading::pool pool(num_workers);
    for (std::size_t ii = 0; ii < num_workers; ++ii)
    {
        pool.emplace(ii, std::ref(targets[ii]));
    }

    auto check_finished = [&]()
    {
        bool all_done = true;
        for (std::size_t jj = 0; jj < num_workers; ++jj)
        {
            all_done = all_done && targets[jj].done;
        }
        return all_done;
    };

    for (std::size_t ii = 0; ii < 10; ++ii)
    {
        pool.notify();
        while (!check_finished())
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        for (std::size_t jj = 0; jj < num_workers; ++jj)
        {
            targets[jj].done = false;
        }
    }

    for (std::size_t ii = 0; ii < iters; ++ii)
    {
        for (std::size_t jj = 0; jj < num_workers; ++jj)
        {
            assert(targets[jj].counter == iters);
        }
    }

    pool.join();
}

int main()
{
    test_worker_notify();
    test_thread_pool();
}
