#include <catch2/catch.hpp>

#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <image/image_gen.hpp>

#include <iostream>

TEST_CASE("Test replicated padding of size 1")
{
    const int colors[2][2 * 3] =
    {
        {
            1, 2, 3,
            4, 5, 6
        },
        {
            6, 5, 4,
            3, 2, 1
        }
    };

    const int expected[2][(2 + 2) * (3 + 2)] =
    {
        {
            1, 1, 2, 3, 3,
            1, 1, 2, 3, 3,
            4, 4, 5, 6, 6,
            4, 4, 5, 6, 6
        },
        {
            6, 6, 5, 4, 4,
            6, 6, 5, 4, 4,
            3, 3, 2, 1, 1,
            3, 3, 2, 1, 1
        }
    };

    sight::Image<short> im(3, 2, 2);
    for (int i = 0, k = 0; i < 2; ++i)
    {
        for (int j = 0; j < 3; ++j, ++k)
        {
            im(i, j, 0) = colors[0][k];
            im(i, j, 1) = colors[1][k];
        }
    }
    
    const int padding = 1;
    const auto res = sight::Pad(im, padding);

    for (int i = 0, k = 0; i < (2 + 2); ++i)
    {
        for (int j = 0; j < (3 + 2); ++j, ++k)
        {
            REQUIRE(res(i, j, 0) == expected[0][k]);
            REQUIRE(res(i, j, 1) == expected[1][k]);
        }
    }
}

TEST_CASE("Test that replicated padding + original roi results in the same image")
{
    const int colors[2][2 * 3] =
    {
        {
            1, 2, 3,
            4, 5, 6
        },
        {
            6, 5, 4,
            3, 2, 1
        }
    };

    sight::Image<short> im(3, 2, 2);
    for (int i = 0, k = 0; i < 2; ++i)
    {
        for (int j = 0; j < 3; ++j, ++k)
        {
            im(i, j, 0) = colors[0][k];
            im(i, j, 1) = colors[1][k];
        }
    }

    // arbitrary padding
    const int padding = 1;
    const auto padded = sight::Pad(im, padding);
    const auto res = padded(sight::Roi(padding, padding, im.w, im.h));

    REQUIRE(res.w == im.w);
    REQUIRE(res.h == im.h);
    REQUIRE(res.c == im.c);

    for (int i = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j)
        {
            for (int c = 0; c < im.c; ++c)
            {
                REQUIRE(res(i, j, c) == im(i, j, c));
            }
        }
    }
}

// TEST_CASE("Test generate corner")
// {
//     // Creates a 2x2 block of squares, where the top-right square
//     // is all white and the rest are black.
//     const int black = 0;
//     const int white = 1;
//     const auto im = sight::CreateCorner<uint8_t>(99, M_PI_4, black, white);
//     const int cx = (im.w + 1) / 2;
//     const int cy = (im.h + 1) / 2;
    
//     for (int i = 0; i < im.h; ++i)
//     {
//         for (int j = 0; j < im.w; ++j)
//         {
//             // In the first quadrant? Should be white
//             if (i < cy && j > cx)
//             {
//                 REQUIRE(im(i, j) == white);
//             }
//             else
//             {
//                 REQUIRE(im(i, j) == black);
//             }
//         }
//     }
// }