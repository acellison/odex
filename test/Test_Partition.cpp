
#include "odex/detail/partition.hpp"
#include <iostream>
#include <vector>

static void print_partitions(std::vector<std::vector<std::size_t>> partitions)
{
    std::cout << "[" << std::endl;
    for (auto const& partition : partitions)
    {
        std::cout << "  [";
        for (std::size_t jj = 0; jj < partition.size(); ++jj)
        {
            std::cout << partition[jj];
            if (jj < partition.size()-1) std::cout << " ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
}

static void test_partition_sized(std::size_t count)
{
    std::vector<std::size_t> indices;;
    for (std::size_t ii = 2; ii <= count; ii += 2)
    {
        indices.push_back(ii);
    }

    auto partitions = odex::detail::partition(indices.begin(), indices.size());
    print_partitions(partitions);
}

static void test_partition()
{
    test_partition_sized(14);
    test_partition_sized(16);

    std::size_t step_counts[] = {2, 16, 18, 20};
    auto partitions = odex::detail::partition(step_counts, 4);
    print_partitions(partitions);
}

int main()
{
    test_partition();
}
