#include <catch2/catch.hpp>

#include <linear/vec.hpp>
#include <geometric/triangulation.hpp>

TEST_CASE("Random triangulation test")
{
    using namespace sight;
    using S = double;
    
    const S theta = 0.1234;
    const SE3<S> T = {
        {
            cos(theta),   0, -sin(theta),
                     0,   1,           0,
            sin(theta),   0,  cos(theta)
        },
        {
            .05, -.47, 0.005
        }
    };

    const Vec3<S> pt0 = {1.07, 0.0521, .875};
    const Vec3<S> pt1 = T * pt0;

    const Vec3<S> r0 = pt0.Normalized();
    const Vec3<S> r1 = pt1.Normalized();

    Vec3<S> estPt;
    REQUIRE(TriangulateL1Angular(r1, r0, T, estPt));

    const S eps = std::numeric_limits<S>::epsilon();
    REQUIRE(estPt(0) == Approx(pt1(0)).margin(eps));
    REQUIRE(estPt(1) == Approx(pt1(1)).margin(eps));
    REQUIRE(estPt(2) == Approx(pt1(2)).margin(eps));
}