#include <cmath>
#include <iostream>

#include <catch2/catch.hpp>
#include <lie/se3.hpp>
#include <lie/lie_ops.hpp>
#include <Eigen/Core>

TEST_CASE("Test SE3 identity")
{
    using namespace sight;
    using S = float;
    const Vec3<S> v = {1, 2, 3};
    const SE3<S> R = SE3<S>::Identity();
    const auto Rv = R * v;

    REQUIRE(Rv(0) == v(0));
    REQUIRE(Rv(1) == v(1));
    REQUIRE(Rv(2) == v(2));
}


TEST_CASE("Test SE3 mult Vec3")
{
    using namespace sight;
    using S = float;

    const SE3<S> T = {
        {
            0, -1, 0,
            1,  0, 0,
            0,  0, 1
        },
        { S(.1), S(-.2), S(.3)}
    };

    const Vec3<S> v = {1, 2, 3};
    const auto Tv = T * v;

    REQUIRE(Tv(0) == Approx(-v(1) + S(.1)));
    REQUIRE(Tv(1) == Approx( v(0) - S(.2)));
    REQUIRE(Tv(2) == Approx( v(2) + S(.3)));
}

TEST_CASE("Test SE3 Inverse")
{
    using namespace sight;
    using S = float;

    const S theta = S(.79812);
    const S cost = cos(theta);
    const S sint = sin(theta);

    const SE3<S> T = {
        {
            cost, -sint, S(0),
            sint,  cost, S(0),
            S(0),  S(0), S(1)
        },
        {1, 2, 3}
    };

    const auto Tinv = T.Inverse();
    const auto I = T * Tinv;

    const double eps = std::numeric_limits<double>::epsilon();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == j)
            {
                REQUIRE(I(i, j) == Approx(1.0).margin(eps));
            }
            else
            {
                REQUIRE(I(i, j) == Approx(0.0).margin(eps));
            }
        }
    }
}

TEST_CASE("Test SE3 Log")
{
    using namespace sight;
    using S = double;

    S rx, ry, rz;

    S ux = S(.1);
    S uy = S(29.7);
    S uz = S(-5.0003);

    SECTION("Test 1")
    {
        // Create random vector of magnitude theta
        const S theta = 1.23456789;
        Vec3<S> v = {S(1), S(2.5), S(3.333)};
        v *= (theta / v.Norm());
        
        rx = v(0);
        ry = v(1);
        rz = v(2);
    }

    SECTION("Test 2")
    {
        // Theta close to 0
        const S theta = S(acos(.9999999));
        Vec3<S> v = {S(1), S(2.5), S(3.333)};
        v = v * (theta / v.Norm());
        
        rx = v(0);
        ry = v(1);
        rz = v(2);
    }

    SECTION("Test 3")
    {
        // Theta close to pi
        const S theta = S(acos(-.9999999));
        Vec3<S> v = {S(1), S(2.5), S(3.333)};
        v *= (theta / v.Norm());
        
        rx = v(0);
        ry = v(1);
        rz = v(2);
    }

    Eigen::Matrix4d A;
    A <<
        0, -rz, ry, ux,
        rz, 0, -rx, uy,
        -ry, rx, 0, uz,
        0, 0, 0, 0;
    
    const Eigen::Matrix4d eA = Exp<S, 4>(A);
    const SE3<S> T = {
        {
            eA(0, 0), eA(0, 1), eA(0, 2),
            eA(1, 0), eA(1, 1), eA(1, 2),
            eA(2, 0), eA(2, 1), eA(2, 2)
        },
        {
            eA(0, 3), eA(1, 3), eA(2, 3)
        }
    };

    const Vec<S, 6> v_est = T.Log();

    const S eps = std::numeric_limits<S>::epsilon();
    REQUIRE(v_est(0) == Approx(rx).margin(eps));
    REQUIRE(v_est(1) == Approx(ry).margin(eps));
    REQUIRE(v_est(2) == Approx(rz).margin(eps));
    REQUIRE(v_est(3) == Approx(ux).margin(eps));
    REQUIRE(v_est(4) == Approx(uy).margin(eps));
    REQUIRE(v_est(5) == Approx(uz).margin(eps));
}
