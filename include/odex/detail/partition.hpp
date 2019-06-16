#ifndef ODEX_DETAIL_PARTITION_HPP
#define ODEX_DETAIL_PARTITION_HPP

#include <algorithm>
#include <iterator>
#include <numeric>
#include <cstddef>
#include <vector>

namespace odex {
namespace detail {

/// Attempt to partition a length n array of integers into k partitions 
/// such that no partition sums to more than maxheight.
/// Assumes a is sorted in descending order
template <class InputIterator>
std::vector<std::vector<std::size_t>> _try_partition(InputIterator a, std::size_t n, std::size_t k, std::size_t maxheight)
{
    using diff_t = typename std::iterator_traits<InputIterator>::difference_type;

    // construct the empty bin vectors
    std::vector<std::vector<std::size_t>> bins(k);

    // current sums of each bin
    std::vector<std::size_t> sums(k, std::size_t(0));

    // used flag for each element in the input data set
    std::vector<bool> used(n);

    // iterate over the bins, placing as much data in each as possible
    for (std::size_t ii = 0; ii < k; ++ii)
    {
        // check all elements
        for (std::size_t jj = 0; jj < n; ++jj)
        {
            // if the element is unused, attempt to place it in the current bin
            if (!used[jj])
            {
                auto tmp = sums[ii]+a[static_cast<diff_t>(jj)];
                if (tmp <= maxheight)
                {
                    // the element fits in the bin.  drop it in and mark it used
                    used[jj] = true;
                    sums[ii] = tmp;
                    bins[ii].push_back(a[static_cast<diff_t>(jj)]);
                }
            }
        }
    }
    
    // check all elements were used
    if (!std::all_of(used.begin(), used.end(), [](bool v){return v;}))
    {
        // not all elements were used.  clear the bins to signal the 
        // partitioning failed
        bins.clear();
    }
    return bins;
}

/// Partition the input data into a vector of bins, each with height
/// no greater than the maximum element in the data.
/// \param a input data to partition
/// \param n number of elements in input data
template <class InputIterator>
std::vector<std::vector<std::size_t>> partition(InputIterator a, std::size_t n)
{
    using diff_t = typename std::iterator_traits<InputIterator>::difference_type;

    // copy the input data and sort in descending order
    using value_type = typename std::iterator_traits<InputIterator>::value_type;
    std::vector<value_type> sorted(a, a+static_cast<diff_t>(n));
    std::sort(sorted.rbegin(), sorted.rend());

    // maxheight is the now first element
    auto maxheight = sorted[0];
    auto sum = std::accumulate(sorted.begin(), sorted.end(), std::size_t(0));

    // ceil divide to find the first number of partitions to check
    auto ceildiv = [](auto num, auto den)
    {
        return std::size_t(-(-int(num)/int(den)));
    };
    auto first = ceildiv(sum, maxheight);

    // loop until we find a partitioning that fits
    for (std::size_t ii = first; ii <= n; ++ii)
    {
        auto bins = _try_partition(sorted.begin(), n, ii, maxheight);
        if (!bins.empty())
        {
            return bins;
        }
    }
    return {};
}

} // namespace odex
} // namespace detail

#endif // ODEX_DETAIL_PARTITION_HPP
