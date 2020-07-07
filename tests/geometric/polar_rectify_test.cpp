// #include <catch2/catch.hpp>

// #include <geometric/polar_rectify.hpp>

// TEST_CASE("Test epipolar region")
// {
//     using namespace sight;

//     Vec2<double> e;
//     int region;
//     SECTION("Region 1")
//     {
//         region = 1;
//         e = {-.5, -.5};
//     }

//     SECTION("Region 2")
//     {
//         region = 2;
//         e = {.25, -.3};
//     }

//     SECTION("Region 3")
//     {
//         region = 3;
//         e = {1.2, -.01};
//     }
    
//     SECTION("Region 4")
//     {
//         region = 4;
//         e = {-.1234, 0.111};
//     }
    
//     SECTION("Region 5")
//     {
//         region = 5;
//         e = {.1, .9};
//     }

//     SECTION("Region 6")
//     {
//         region = 6;
//         e = {1.0001, .9999};
//     }

//     SECTION("Region 7")
//     {
//         region = 7;
//         e = {-10.0, 1.001};
//     }

//     SECTION("Region 8")
//     {
//         region = 8;
//         e = {0.2, 2.5};
//     }
    
//     SECTION("Region 9")
//     {
//         region = 9;
//         e = {2.0, 2.0};
//     }

//     REQUIRE(ComputeEpipolarRegion(e) == region);
// }

// TEST_CASE("Test extermal corners")
// {
//     using namespace sight;

//     Vec2<double> e;
//     Vec2<double> c0, c1;

//     SECTION("Inside image")
//     {
//         e = {.5, .5};
//         REQUIRE(GetExtremalCorners(ComputeEpipolarRegion(e), c0, c1) == false);
//     }

//     SECTION("Left-center")
//     {
//         e = {-.5, .5};
//         REQUIRE(GetExtremalCorners(ComputeEpipolarRegion(e), c0, c1) == true);
//         REQUIRE(c0(0) == 0.0);
//         REQUIRE(c0(1) == 0.0);
//         REQUIRE(c1(0) == 0.0);
//         REQUIRE(c1(1) == 1.0);
//     }
// }

// TEST_CASE("Test min/max angles")
// {
//     using namespace sight;

//     Vec2<double> e;
//     double t0, t1;

//     const double eps = std::numeric_limits<double>::epsilon();

//     SECTION("Inside image")
//     {
//         e = {.5, .5};
//         GetPolarMinMaxAngles(e, t0, t1);
//         REQUIRE(t0 == 0.0);
//         REQUIRE(t1 == 2 * M_PI);
//     }

//     SECTION("Left-center")
//     {
//         e = {-.5, .5};
//         GetPolarMinMaxAngles(e, t0, t1);
//         REQUIRE(t0 == Approx(-M_PI_4).margin(eps));
//         REQUIRE(t1 == Approx(M_PI_4).margin(eps));
//     }

//     SECTION("Bottom-right")
//     {
//         e = {2.0, 1.0};
//         GetPolarMinMaxAngles(e, t0, t1);
//         REQUIRE(t0 == Approx(M_PI).margin(eps));
//         REQUIRE(t1 == Approx(M_PI * 1.25).margin(eps));
//     }
// }

// TEST_CASE("Test min/max rho")
// {
//     using namespace sight;

//     using S = float;
    
//     Vec3<S> _;

//     Vec2<S> e;
//     S minR, maxR;

//     Vec3<S> borders[4];
//     GetImageBorderLines(borders);

//     const S eps = std::numeric_limits<S>::epsilon();
    
//     SECTION("Horizontal mid-line")
//     {
//         e = {S(-.5), S(.5)};
//         GetPolarMinMaxRho(e, { S(1), S(0) }, borders, _, minR, maxR);

//         REQUIRE(minR == S(.5));
//         REQUIRE(maxR == S(1.5));
//     }    

//     SECTION("Point inside image to corner")
//     {
//         e = {S(.25), S(.25)};

//         GetPolarMinMaxRho(e, { S(.5), S(.5) }, borders, _, minR, maxR);

//         REQUIRE(minR == S(0));
//         REQUIRE(maxR == Approx(1.0606601717798212).margin(eps));
//     }

//     SECTION("Point out to infinity")
//     {
//         e = {S(-.5), S(-.5)};
//         GetPolarMinMaxRho(e, { S(-.5), S(-.5) }, borders, _, minR, maxR);
//         REQUIRE(minR == std::numeric_limits<S>::max());
//         REQUIRE(maxR == std::numeric_limits<S>::min());
//     }
// }

//     // Returns an integer in the range [1, 9] representing
//     // where the epipole lies around [0, 1] x [0, 1].
//     //
//     //    1 | 2 | 3
//     //  ----a---b----
//     //    4 | 5 | 6
//     //  ----c---d----
//     //    7 | 8 | 9

//     //   7 | 8 | 9
//     //  ---c---d---
//     //   4 | 5 | 6
//     //  ---a---b---
//     //   1 | 2 | 3