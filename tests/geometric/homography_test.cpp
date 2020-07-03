#include <catch2/catch.hpp>

#include <cmath>
#include <iterator>
#include <geometric/homography_estimation.hpp>

TEST_CASE("Make 2d vectors into 3d homogeneous vectors")
{
    Eigen::Vector2d u(10.1234, -573.2);
    auto v = sight::ToHomogeneous(u);

    REQUIRE(v(0) == u(0));
    REQUIRE(v(1) == u(1));
    REQUIRE(v(2) == 1.0);
}

TEST_CASE("Make 3d vectors into 2d homogeneous vectors")
{
    const double x = -5.4321;
    const double y = 2975.0;
    const double s = 1.234;

    Eigen::Vector3d u(s * x, s * y, s);
    auto v = sight::FromHomogeneous(u);

    REQUIRE(v(0) == Approx(x));
    REQUIRE(v(1) == Approx(y));
}

TEST_CASE("Apply scale + translational homography")
{
    const double s1 = 2.2;
    const double s2 = cos(.5);
    const double tx = .12345;
    const double ty = -5.4321;
    Eigen::Matrix3d H;
    H <<
        s1, 0, tx,
        0, s2, ty,
        0, 0, 1;
    
    const Eigen::Vector2d xy(10.5, 20.579);
    const auto v = sight::ApplyHomography(H, xy);

    REQUIRE(v(0) == Approx(s1 * xy(0) + tx));
    REQUIRE(v(1) == Approx(s2 * xy(1) + ty));
}

TEST_CASE("Normalization centers points to origin and scales to sqrt(2)")
{
    // Create a quad of points
    sight::EigenList<Eigen::Vector2d> vs;
    vs.push_back({1.2, -0.7});
    vs.push_back({10.791, 5.3});
    vs.push_back({2.2, 2.2});
    vs.push_back({6, 0});

    // Find normalization matrix
    auto T = sight::NormalizationMatrix(vs);

    // Compute mean and average magnitude
    Eigen::Vector2d mean = Eigen::Vector2d::Zero();
    double sum_mag = 0.0;
    for (const auto& v : vs)
    {
        const auto Tv = T * Eigen::Vector3d(v(0), v(1), 1);
        mean(0) += Tv(0);
        mean(1) += Tv(1);
        sum_mag += sqrt(Tv(0) * Tv(0) + Tv(1) * Tv(1));
    }
    mean *= (1.0 / vs.size());
    sum_mag /= vs.size();

    // Check conditions
    const double eps = std::numeric_limits<double>::epsilon();
    REQUIRE(mean(0) == Approx(0.0).margin(eps));
    REQUIRE(mean(1) == Approx(0.0).margin(eps));
    REQUIRE(sum_mag == Approx(sqrt(2)));
}

TEST_CASE("Inversion of scaling + translation matrix works")
{
    // Create an arbitary scale + translation matrix
    const double s = 10.32;
    const double x = -.054321;
    const double y = 31.7978;
    Eigen::Matrix3d T;
    T <<
        s, 0, s * x,
        0, s, s * y,
        0, 0, 1;
    
    // Invert it quickly
    const auto Tinv = sight::InvertScaleAndTrans(T);

    // The expected result of multiplying matrices which
    // are inverses of one another is the identity matrix.
    const auto expected = Eigen::Matrix3d::Identity();
    const auto actual = T * Tinv;

    const double eps = std::numeric_limits<double>::epsilon();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            REQUIRE(actual(i, j) == Approx(expected(i, j)).margin(eps));
        }
    }
}

#include <iostream>
TEST_CASE("Homography maps similarity transform")
{
    using namespace std;
    const auto s = 1.23456789;
    const auto th = .2;
    const auto tx = -30.5;
    const auto ty = 100.3;
    Eigen::Matrix3d T;
    T <<
        s * cos(th),    -sin(th), 0,
            sin(th), s * cos(th), 0,
                  0,           0, 1;

    sight::EigenList<Eigen::Vector2d> src = {
        {0, 0},
        {0, 1},
        {1, 1},
        {1, 0}
    };

    sight::EigenList<Eigen::Vector2d> dst(src.size());
    for (int i = 0; i < dst.size(); ++i)
    {
        dst[i] = sight::ApplyHomography(T, src[i]);
    }

    Eigen::Matrix3d H;
    SECTION("with point normalization")
    {
        H = sight::FitHomography(src, dst, true);
    }
    
    SECTION("without point normalization")
    {
        H = sight::FitHomography(src, dst, false);
    }

    const auto eps = std::numeric_limits<double>::epsilon();
    for (int i = 0; i < src.size(); ++i)
    {
        const auto& expected = dst[i];
        const auto actual = sight::ApplyHomography(H, src[i]);

        REQUIRE(actual(0) == Approx(expected(0)).margin(eps));
        REQUIRE(actual(1) == Approx(expected(1)).margin(eps));
    }
}

TEST_CASE("Ransac correctly guesses the inlier set")
{
    // Create a random homography
    Eigen::Matrix3d H;
    H <<
        0.57, .750,  .001,
        -3.9,  .57,  -.28,
        .866, .538, .1235;

    sight::EigenList<Eigen::Vector2d> src, dst;

    const double eps = 1e-3;
    const double max_noise = eps / 10;
    std::default_random_engine gen;
    std::normal_distribution<double> dist(0.0, eps * eps);

    const int rows = 12;
    const int cols = 12;
    const int nGood = 2 * (rows * cols) / 3;
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            src.push_back({i, j});

            if ((i + j) % 3 != 0)
            {
                // Make a good mapping with a tiny bit of
                // capped gaussian noise.
                dst.push_back(sight::ApplyHomography(H, src.back()));
                dst.back()(0) += std::clamp(-max_noise, dist(gen), max_noise);
                dst.back()(1) += std::clamp(-max_noise, dist(gen), max_noise);
            }
            else
            {
                // Make a bad mapping (H must be I, which it's not)
                dst.push_back({i, j});
            }
        }
    }

    std::vector<int> inliers;
    const auto H_est = sight::RansacHomography(src, dst, 100, eps, inliers);

    REQUIRE(inliers.size() >= nGood * .9);

    for (int i : inliers)
    {
        const auto& expected = dst[i];
        const auto actual = sight::ApplyHomography(H_est, src[i]);
        
        REQUIRE(actual(0) == Approx(expected(0)).margin(eps));
        REQUIRE(actual(1) == Approx(expected(1)).margin(eps));
    }
}