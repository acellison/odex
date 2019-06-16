
#ifndef ODEX_CONVECTOR_HPP
#define ODEX_CONVECTOR_HPP

#include "central_difference.hpp"

template <class T, class Matrix>
class convector
{
public:
    using value_type = T;
    using matrix_type = Matrix;

    convector(value_type k, value_type cx, value_type cy)
    : m_k(k), m_cx(cx), m_cy(cy)
    {    }

    auto operator()(value_type, matrix_type const& u)
    {
        central_difference(u, m_k, m_ux, m_uy);
        return m_cx*m_ux+m_cy*m_uy;
    }

private:
    matrix_type m_ux;
    matrix_type m_uy;
    value_type  m_k;
    value_type  m_cx;
    value_type  m_cy;
};

#endif // ODEX_CONVECTOR_HPP
