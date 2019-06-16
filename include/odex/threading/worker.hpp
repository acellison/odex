
#ifndef ODEX_THREADING_WORKER_HPP
#define ODEX_THREADING_WORKER_HPP

#include "odex/threading/semaphore.hpp"
#include <utility>
#include <atomic>
#include <thread>

namespace odex {
namespace threading {

/// Worker thread that waits for notification before processing data.  The
/// worker has its target load attached at construction, and will call this
/// each time notify() is called, then go back to waiting on a semaphore.
/// When passing function and arguments to the worker thread, copies are made
/// to avoid storing references to stack data that gets destroyed after function
/// exit.  To pass function object and arguments via reference, wrap the 
/// arguments to the constructor with a std::reference_wrapper via std::ref().
class worker
{
    worker(worker const&) = delete;
public:
    /// Construct the worker with a function and arguments to pass on to it.
    /// Arguments are passed on to the worker thread by value to prevent bad
    /// access on objects allocated on the caller's stack.  Use std::ref() to
    /// pass by reference if no copies are desired.
    template <class Function, class... Args, 
              class = std::enable_if_t<!std::is_same<std::decay_t<Function>, worker>{}>>
    explicit worker(Function&& function, Args&&... args)
    : m_exit_flag(false)
    , m_semaphore()
    , m_thread(_make_target(std::forward<Function>(function), std::forward<Args>(args)...))
    {

    }

    /// Destroy the worker, joining.
    ~worker()
    {
        join();
    }

    /// Notify the worker that data is ready to be processed.
    void notify()
    {
        m_semaphore.notify();
    }

    /// Join the worker before destruction.
    void join()
    {
        if (m_thread.joinable())
        {
            m_exit_flag = true;
            notify();
            m_thread.join();
        }
    }

private:
    /// Make the target function.
    template <class Function, class... Args>
    auto _make_target(Function&& function, Args&&... args)
    {
        // Captures function and arguments by value
        return [this, function, args...]()
        {
            _run(function, args...);
        };
    }

    /// Main run loop.
    template <class Function, class... Args>
    void _run(Function&& function, Args&&... args)
    {
        while (true)
        {
            // Wait for notification
            m_semaphore.wait();

            // If exit signal, break from the loop
            if (m_exit_flag) break;

            // Call the target function
            std::forward<Function>(function)(std::forward<Args>(args)...);
        }
    }

private:
    std::atomic<bool> m_exit_flag;
    semaphore m_semaphore;
    std::thread m_thread;
};

} // namespace threading
} // namespace odex

#endif // ODEX_THREADING_WORKER_HPP
