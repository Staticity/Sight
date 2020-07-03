#include <limits>
#include <iostream>

#include <catch2/catch.hpp>
#include <Eigen/Core>
#include <lie/exp_utils.hpp>

TEST_CASE("Test Computation of Exp Matrix 3x3")
{
    // Should compute I + a[w]x + b[w]x^2

    using S = float;
    S x, y, z;
    S a, b;

    SECTION("Test 1")
    {
        x = S(.1234);
        y = S(.0890637);
        z = S(-5.3);

        a = cos(S(2.7));
        b = S(-.9);
    }

    Eigen::Matrix3<double> wx;
    wx <<
         0, -z,  y,
         z,  0, -x,
        -y,  x,  0;

    const Eigen::Matrix3<double> actual =
        Eigen::Matrix3<double>::Identity() + a * wx + b * (wx * wx);

    S M[3 * 3];
    sight::ComputeExpMatrix3x3(&M[0], a, b, x, y, z);

    const S eps = std::numeric_limits<S>::epsilon();
    for (int i = 0, k = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j, ++k)
        {
            REQUIRE(M[k] == Approx(actual(i, j)).margin(eps));
        }
    }
}