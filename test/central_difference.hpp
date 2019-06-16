
#ifndef ODEX_CENTRAL_DIFFERENCE_HPP
#define ODEX_CENTRAL_DIFFERENCE_HPP

#include <cstddef>

template <class Matrix, class GridSpacing>
auto central_difference(Matrix const& u, GridSpacing k, Matrix& ux, Matrix& uy)
{
    std::ptrdiff_t n = u.rows();
    std::ptrdiff_t m = u.cols();
    for (std::ptrdiff_t ii = 0; ii < n; ++ii)
    {
        for (std::ptrdiff_t jj = 1; jj < m-1; ++jj)
        {
            ux(ii,jj) = (u(ii,jj+1)-u(ii,jj-1))/(2*k);
        }
        ux(ii,0)   = (u(ii,1)-u(ii,m-1))/(2*k);
        ux(ii,m-1) = (u(ii,0)-u(ii,m-2))/(2*k);
    }
    for (std::ptrdiff_t ii = 0; ii < m; ++ii)
    {
        for (std::ptrdiff_t jj = 1; jj < n-1; ++jj)
        {
            uy(ii,jj) = (u(jj+1,ii)-u(jj-1,ii))/(2*k);
        }
        uy(0,ii)   = (u(1,ii)-u(n-1,ii))/(2*k);
        uy(n-1,ii) = (u(0,ii)-u(n-2,ii))/(2*k);
    }
}

#endif // ODEX_CENTRAL_DIFFERENCE_HPP
