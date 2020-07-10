#include "camera_test_helpers.hpp"

#include <catch2/catch.hpp>

#include <calibration/equidistantmodel.hpp>
#include <calibration/pinholemodel.hpp>
#include <calibration/radial4model.hpp>
#include <calibration/radialtanmodel.hpp>
#include <linear/vec.hpp>

TEST_CASE("Test Pinhole")
{
    using namespace sight;
    using S = double;
    const PinholeModel<S> cam(
        S(.7), S(.4222), S(.5703), S(.45123));

    TestCamera(cam);
}

TEST_CASE("Test Radial4")
{
    using namespace sight;
    using S = double;
    const Radial4Model<S> cam(
        S(.8), S(.4), S(.512), S(.4933), S(.1), S(.02), S(-.03), S(.0045));

    TestCamera(cam);
}

TEST_CASE("Test Equidistant")
{
    using namespace sight;
    using S = double;
    const EquidistantModel<S> cam(
        S(380.81042871360756),
        S(380.81194179427075),
        S(510.29465304840727),
        S(514.33046305385062),
        S(0.010171079892421483),
        S(-0.010816440029919381),
        S(0.0059427817694127560),
        S(-0.0016622846678576431));

    TestCamera(cam);
}

TEST_CASE("Test RadialTan")
{
    using namespace sight;
    using S = double;
    const RadialTanModel<S> cam(
        S(935.1344941116828),
        S(931.4115487638016),
        S(639.5005137221706),
        S(338.0236559548098),
        S(-.001),
        S(.001),
        S(.002),
        S(.003),
        S(.004),
        S(.005),
        S(.006),
        S(.007));

    TestCamera(cam);
}
