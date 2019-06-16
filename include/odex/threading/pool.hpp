
#ifndef ODEX_THREADING_POOL_HPP
#define ODEX_THREADING_POOL_HPP

#include "odex/threading/worker.hpp"
#include "odex/threading/semaphore.hpp"
#include <cassert>
#include <cstddef>
#include <vector>
#include <memory>

namespace odex { 
namespace threading {

/// Thread pool object that manages a number of worker threads.  Each worker
/// may have its own distinct target load function.  Individual workers can
/// be told to process via calls to notify().  To guarantee synchronization,
/// process() dispatches all the workers and does not return until all have
/// completed their work.
class pool
{
    pool(pool const&) = delete;
public:
    /// Construct the thread pool with a number of workers
    explicit pool(std::size_t num_workers)
    : m_workers(num_workers)
    , m_completion_count(0)
    {

    }

    /// Destroy the thread pool, joining the worker threads
    ~pool()
    {
        join();
    }

    /// Number of workers in the pool.
    std::size_t size() const
    {
        return m_workers.size();
    }

    /// Construct a worker in place.  Function and args are passed by value to
    /// the worker thread to prevent usage of stack-allocated objects after
    /// their destruction.  To pass by reference, wrap the function object
    /// and arguments in a std::reference_wrapper via std::ref().
    template <class Function, class... Args>
    void emplace(std::size_t index, Function&& function, Args&&... args)
    {
        assert(index < m_workers.size() && "Worker index out of range!");
        m_workers[index].reset(new worker(_make_target(std::forward<Function>(function), std::forward<Args>(args)...)));
    }

    /// Tell the workers to process, synchronizing.  This call does not return
    /// until all workers have finished processing.  If other work can be done
    /// concurrently on the calling thread, use notify() instead.  Managing
    /// worker processing completion is then left to the caller.
    void process()
    {
        std::size_t count = m_workers.size();
        for (std::size_t ii = 0; ii < count; ++ii)
        {
            notify(ii);
        }
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cv.wait(lock, [this, count]{return m_completion_count == count;});
        m_completion_count = 0;
    }

    /// Notify all worker to process.
    void notify()
    {
        std::size_t count = m_workers.size();
        for (std::size_t ii = 0; ii < count; ++ii)
        {
            notify(ii);
        }
    }

    /// Notify the worker at index to process.
    void notify(std::size_t index)
    {
        if (m_workers[index])
        {
            m_workers[index]->notify();
        }
    }

    /// Join all workers to this thread.
    void join()
    {
        std::size_t count = m_workers.size();
        for (std::size_t ii = 0; ii < count; ++ii)
        {
            join(ii);
        }
    }

    /// Join the worker at index to this thread.
    void join(std::size_t index)
    {
        if (m_workers[index])
        {
            m_workers[index]->join();
        }
    }

private:
    /// Make the target function.
    template <class Function, class... Args>
    auto _make_target(Function&& function, Args&&... args)
    {
        return [this, function, args...]()
        {
            // Call the target function
            function(args...);

            // Get the lock
            std::lock_guard<std::mutex> lock(m_mutex);

            // Bump the number of finished threads
            ++m_completion_count;

            // If all workers are done, notify the main thread
            if (m_completion_count == m_workers.size())
            {
                m_cv.notify_one();
            }
        };
    }

private:
    std::vector<std::unique_ptr<worker>> m_workers;
    std::condition_variable m_cv;
    std::mutex m_mutex;
    std::size_t m_completion_count;
};

} // namespace odex
} // namespace threading

#endif // ODEX_THREADING_POOL_HPP
