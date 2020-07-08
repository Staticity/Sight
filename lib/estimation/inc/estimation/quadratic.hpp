#pragma once

#include <Eigen/LU>

namespace sight
{
    template <typename S>
    class Quadratic
    {
    public:

        Quadratic()
        {}

        Quadratic(
            S x0, S y0,
            S x1, S y1,
            S x2, S y2)
        {
            // Find ax{i}^2 + bx{i} + c = y{i}
            //
            // for all (x{i}, y{i}).

            // [x0^2  x0  1] [a]   [y0]
            // [x1^2  x1  1] [b] = [y1]
            // [x2^2  x2  1] [c]   [y2]
            //
            // Perform least squares estimate with pseudo-inverse
            Fit(x0, y0, x1, y1, x2, y2, a, b, c);
        }

        S Eval(const S x) const
        {
            return a * x * x + b * x + c;
        }

        bool GetExtremum(S& x) const
        {
            // Find where the derivative of the expression is equal to 0.
            //
            //     d/dx(ax^2 + bx + c) = 0
            //
            //     2ax + b = 0
            if (abs(a) < std::numeric_limits<S>::epsilon())
            {
                return false;
            }
            
            x = -b / (2 * a);
            return true;
        }

        bool HasMaximum() const
        {
            return a > std::numeric_limits<S>::epsilon();
        }

        bool HasMinimum() const
        {
            return a < -std::numeric_limits<S>::epsilon();
        }

        static void Fit(
            const S x0, const S y0,
            const S x1, const S y1,
            const S x2, const S y2,
            S& a, S& b, S& c)
        {
            Eigen::Matrix3<S> A;
            A <<
                x0 * x0, x0, 1,
                x1 * x1, x1, 1,
                x2 * x2, x2, 1;
            
            Eigen::Vector3<S> v(y0, y1, y2);
            Eigen::Vector3<S> x = (A.transpose() * A).inverse() * (A.transpose() * v);
            a = x(0);
            b = x(1);
            c = x(2);
        }

        static void Fit(
            const S* xs,
            const S* ys,
            const int n,
            S& a, S& b, S& c)
        {
            Eigen::MatrixX<S> A;
            Eigen::VectorX<S> v;
            A.resize(n, 3);
            v.resize(n);

            for (int i = 0; i < n; ++i)
            {
                A.row(i) << xs[i] * xs[i], xs[i], S(1);
                v(i) = ys[i];
            }

            const Eigen::Vector3<S> x = (A.transpose() * A).inverse() * (A.transpose() * v);
            a = x(0);
            b = x(1);
            c = x(2);
        }

        S a, b, c;
    };
}