#pragma once

#include <Eigen/Core>
#include <limits>

namespace sight
{
    template <typename S, int N>
    Eigen::Matrix<S, N, N> Exp(
        const Eigen::Matrix<S, N, N>& A,
        const int iterations = -1,
        const S eps = std::numeric_limits<double>::epsilon())
    {
        // A^k
        Eigen::Matrix<S, N, N> Ak = Eigen::Matrix<S, N, N>::Identity();

        // X = I + A + A^2/2! + A^3/3!
        auto X = Ak; // X = I to start.

        int it = 1;
        S scale = S(1);
        while (iterations == -1 || it <= iterations)
        {
            Ak *= A;
            const Eigen::Matrix<S, N, N> B = Ak * scale;

            if (B.squaredNorm() < eps)
            {
                break;
            }

            X += B;
            
            // The iteration counter manages the factorial
            // factor as well.
            ++it;
            scale /= it;
        }

        return X;
    }
}