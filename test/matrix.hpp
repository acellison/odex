
#ifndef ODEX_MATRIX_HPP
#define ODEX_MATRIX_HPP

#include <Eigen/Core>

/// Dynamically allocated Eigen::Matrix object with default constructor that sets size
/// appropriately.  Eigen static_asserts if large matrices are instantiated on the stack,
/// but we avoid plumbing the resize() method through odex by creating default-constructible
/// matrices that do the work for us.
template <class T, int Rows, int Cols>
class matrix : public Eigen::Matrix<T,-1,-1>
{
    using base = Eigen::Matrix<T,-1,-1>;
public:
    // This constructor allows you to construct matrix from Eigen expressions
    template <class OtherDerived>
    matrix(const Eigen::MatrixBase<OtherDerived>& other)
    : base(other)
    {    }

    // This method allows you to assign Eigen expressions to matrix
    template <class OtherDerived>
    matrix& operator=(Eigen::MatrixBase <OtherDerived> const& other)
    {
        this->base::operator=(other);
        return *this;
    }

    matrix()
    : base(Rows, Cols)
    {    }
};

#endif // ODEX_MATRIX_HPP
