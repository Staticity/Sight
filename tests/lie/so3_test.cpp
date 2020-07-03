#include <cmath>
#include <iostream>

#include <catch2/catch.hpp>
#include <lie/so3.hpp>
#include <lie/lie_ops.hpp>
#include <Eigen/Core>

TEST_CASE("Test SO3 identity")
{
    using namespace sight;
    using S = float;
    const Vec3<S> v = {1, 2, 3};
    const SO3<S> R = SO3<S>::Identity();
    const auto Rv = R * v;

    REQUIRE(Rv(0) == v(0));
    REQUIRE(Rv(1) == v(1));
    REQUIRE(Rv(2) == v(2));
}


TEST_CASE("Test SO3 mult Vec3")
{
    using namespace sight;
    using S = float;

    const SO3<S> R = {
        0, -1, 0,
        1,  0, 0,
        0,  0, 1
    };

    const Vec3<S> v = {1, 2, 3};
    const auto Rv = R * v;

    REQUIRE(Rv(0) == Approx(-v(1)));
    REQUIRE(Rv(1) == Approx(v(0)));
    REQUIRE(Rv(2) == Approx(v(2)));
}

TEST_CASE("Test SO3 Inverse")
{
    using namespace sight;
    using S = float;

    const S theta = S(.79812);
    const S cost = cos(theta);
    const S sint = sin(theta);

    const SO3<S> R = {
        cost, -sint, S(0),
        sint,  cost, S(0),
        S(0),  S(0), S(1)
    };

    const auto Rinv = R.Inverse();
    const auto I = R * Rinv;

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

TEST_CASE("Test SO3 Exponentiation")
{
    using S = double;
    const S x =  1.5 * S(3.141592653);
    const S y =  S(.000023);
    const S z = S(-.000045);

    Eigen::Matrix3<S> A;
    A <<
         0, -z,  y,
         z,  0, -x,
        -y,  x,  0;
    
    const Eigen::Matrix3<S> expA = sight::Exp<S, 3>(A);
    const sight::SO3<S> R = sight::SO3<S>::Exp({x, y, z});

    const S eps = sqrt(std::numeric_limits<S>::epsilon());
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            REQUIRE(R(i, j) == Approx(expA(i, j)).margin(eps));
        }
    }
}

TEST_CASE("Test SO3 Log")
{
    using namespace sight;
    using S = double;

    S x, y, z;

    SECTION("Test 1")
    {
        // Create random vector of magnitude theta
        const S theta = 1.23456789;
        Vec3<S> v = {S(1), S(2.5), S(3.333)};
        v *= (theta / v.Norm());
        
        x = v(0);
        y = v(1);
        z = v(2);
    }

    SECTION("Test 2")
    {
        // Theta close to 0
        const S theta = S(acos(.9999999));
        Vec3<S> v = {S(1), S(2.5), S(3.333)};
        v = v * (theta / v.Norm());
        
        x = v(0);
        y = v(1);
        z = v(2);
    }

    SECTION("Test 3")
    {
        // Theta close to pi
        const S theta = S(acos(-.9999999));
        Vec3<S> v = {S(1), S(2.5), S(3.333)};
        v *= (theta / v.Norm());
        
        x = v(0);
        y = v(1);
        z = v(2);
    }

    Eigen::Matrix3d A;
    A <<
        0, -z, y,
        z, 0, -x,
        -y, x, 0;
    
    const Eigen::Matrix3d eA = Exp<S, 3>(A);
    const SO3<S> R = {
        eA(0, 0), eA(0, 1), eA(0, 2),
        eA(1, 0), eA(1, 1), eA(1, 2),
        eA(2, 0), eA(2, 1), eA(2, 2)
    };

    const Vec3<S> v_est = R.Log();

    const S eps = std::numeric_limits<S>::epsilon();
    REQUIRE(v_est(0) == Approx(x).margin(eps));
    REQUIRE(v_est(1) == Approx(y).margin(eps));
    REQUIRE(v_est(2) == Approx(z).margin(eps));
}
