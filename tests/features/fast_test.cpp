#include <catch2/catch.hpp>

#include <iostream>
#include <image/image.hpp>
#include <image/image_gen.hpp>
#include <features/fast.hpp>

// TEST_CASE("Test fast feature angle calculation")
// {
//     // This generates a 2x2 block image where the upper-left
//     // block is white and the others are black. The angle generated
//     // should be pointing in the (+1, +1) quadrant producing an
//     // angle .25 * M_PI.
//     double expectedAngle;
    
//     SECTION("Angle of 45 degrees")
//     {
//         expectedAngle = M_PI_4;
//     }


//     SECTION("Angle of 30 degrees")
//     {
//         expectedAngle = M_PI / 6;
//     }

//     const auto im = sight::CreateCorner<uint8_t>(99, expectedAngle, 0, 1);
//     const float angle = sight::AngleFAST(im, 50.f, 50.f);
//     REQUIRE(angle == Approx(expectedAngle));
// }
