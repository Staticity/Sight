#include <catch2/catch.hpp>

#include <linear/vec.hpp>
#include <geometric/triangulation.hpp>

TEST_CASE("Random triangulation test")
{
    using namespace sight;
    using S = double;
    const S theta = 0.1234;
    const Vec3<S> pt = {1.07, 0.0521, -.875};
    const SE3<S> T = {
        {
            cos(theta),   0, -sin(theta),
                     0,   1,           0,
            sin(theta),   0,  cos(theta)
        },
        {
            1.234, -.47, 0.0005
        }
    };

    const Vec3<S> r0 = pt.Normalized();
    const Vec3<S> r1 = (T * pt).Normalized();

    Vec3<S> estPt;
    REQUIRE(TriangulateL1Angular(r0, r1, T, estPt));

    const S eps = std::numeric_limits<S>::epsilon();
    REQUIRE(estPt(0) == Approx(pt(0)).margin(eps));
    REQUIRE(estPt(1) == Approx(pt(1)).margin(eps));
    REQUIRE(estPt(2) == Approx(pt(2)).margin(eps));
}