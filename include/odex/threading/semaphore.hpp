
#ifndef ODEX_THREADING_SEMAPHORE_HPP
#define ODEX_THREADING_SEMAPHORE_HPP

#include <condition_variable>
#include <mutex>

namespace odex {
namespace threading {

/// Simple semaphore class for synchronizing worker threads.  One thread waits
/// for notification by calling the wait() method, while another thread can 
/// notify() the waiter that data is available to be processed.
class semaphore
{
public:
    /// Construct the semaphore.
    semaphore()
    : m_ready(false)
    {    }
    
    /// Notify one thread waiting on the semaphore.
    void notify()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_ready = true;
        lock.unlock();
        m_cv.notify_one();
    }
    
    /// Wait for notification.
    void wait()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cv.wait(lock, [this]{return m_ready;});
        m_ready = false;
        lock.unlock();
    }

private:
    std::condition_variable m_cv;
    std::mutex m_mutex;
    bool m_ready;
};

} // end namespace threading
} // end namespace odex

#endif // ODEX_THREADING_SEMAPHORE_HPP

