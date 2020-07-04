#include <catch2/catch.hpp>

#include <iostream>
#include <numeric>

#include <image/image_ops.hpp>
#include <image/image_gen.hpp>

TEST_CASE("Test Gaussian Kernel weights sum to 1")
{
    const std::vector<int> radii = {1, 2, 5, 10};
    const std::vector<float> sigmas = {.25, .5, 1, 3};
    for (int r : radii)
    {
        for (float s : sigmas)
        {
            const auto& kern = sight::GaussianKernel(r, s);
            const auto sum = std::accumulate(kern.begin(), kern.end(), 0.f);

            REQUIRE(sum == Approx(1.0));
        }
    }
}

TEST_CASE("Test 1D Horizontal convolution")
{
    using namespace sight;

    const int colors[2 * 3][3] =
    {
        {1, 2, 3}, {4, 5, 6}, {7, 8, 9},
        {5, 5, 5}, {6, 7, 8}, {5, 5, 5}
    };

    const int expected_colors[2 * 3][3] =
    {
        {3, 3, 3}, {6, 6, 6}, { 3,  3,  3},
        {1, 2, 3}, {0, 0, 0}, {-1, -2, -3}
    };

    Image<int> im(3, 2, 3);
    for (int i = 0, k = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j, ++k)
        {
            im(i, j, 0) = colors[k][0];
            im(i, j, 1) = colors[k][1];
            im(i, j, 2) = colors[k][2];
        }
    }

    const std::vector<short> kernel = {-1, 0, 1};

    const int p = int(kernel.size() / 2);
    const auto view = Pad(im, p)(Roi(p, p, im.w, im.h));
    const auto res = ConvolveHorizontal1D<int, short, int>(view, kernel);
    
    for (int i = 0, k = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j, ++k)
        {
            REQUIRE(res(i, j, 0) == expected_colors[k][0]);
            REQUIRE(res(i, j, 1) == expected_colors[k][1]);
            REQUIRE(res(i, j, 2) == expected_colors[k][2]);
        }
    }
}

TEST_CASE("Test 1D Vertical convolution")
{
    using namespace sight;

    const int colors[3 * 2][3] =
    {
        {1, 2, 3}, {5, 5, 5},
        {4, 5, 6}, {6, 7, 8},
        {7, 8, 9}, {5, 5, 5}
    };

    const int expected_colors[3 * 2][3] =
    {
        {3, 3, 3}, { 1,  2,  3},
        {6, 6, 6}, { 0,  0,  0},
        {3, 3, 3}, {-1, -2, -3}
    };

    Image<int> im(2, 3, 3);
    for (int i = 0, k = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j, ++k)
        {
            im(i, j, 0) = colors[k][0];
            im(i, j, 1) = colors[k][1];
            im(i, j, 2) = colors[k][2];
        }
    }

    const std::vector<short> kernel = {-1, 0, 1};
    const int p = int(kernel.size() / 2);
    const auto view = Pad(im, p)(Roi(p, p, im.w, im.h));
    const auto res = ConvolveVertical1D<int, short, int>(view, kernel);
    
    for (int i = 0, k = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j, ++k)
        {
            REQUIRE(res(i, j, 0) == expected_colors[k][0]);
            REQUIRE(res(i, j, 1) == expected_colors[k][1]);
            REQUIRE(res(i, j, 2) == expected_colors[k][2]);
        }
    }
}

TEST_CASE("Test 2D Convolution with Separable Filter")
{
    // Our 2D Kernel:
    //
    // [1  4  6  4 1]
    // [4 16 24 16 4]
    // [6 24 36 24 6]
    // [4 16 24 16 4]
    // [1  4  6  4 1]
    //
    // Can be produced by:
    //
    //      v^T * v
    //
    // where v = [1 4 6 4 1]
    const std::vector<int> kernel = {1, 4, 6, 4, 1};
    
    const int colors[3 * 3][3] = {
        {1, 2, 3}, {1, 2, 3}, {1, 2, 3},
        {1, 2, 3}, {1, 2, 3}, {1, 2, 3},
        {1, 2, 3}, {1, 2, 3}, {1, 2, 3}
    };
    
    sight::Image<int> im(3, 3, 3);
    for (int i = 0, k = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j, ++k)
        {
            for (int ch = 0; ch < im.c; ++ch)
            {
                im(i, j, ch) = colors[k][ch];
            }
        }
    }

    const int p = int(kernel.size() / 2);
    const auto padded = Pad(im, p);
    const auto view = padded(sight::Roi(p, p, im.w, im.h));
    const auto res = sight::ConvolveSeparable<int, int, int>(view, kernel, kernel);

    // The first channel is all 1's, so the result is simply the sum of the 2D kernel weights
    const int sumWeights = 256;
    for (int i = 0, k = 0; i < im.h; ++i)
    {
        for (int j = 0; j < im.w; ++j, ++k)
        {
            for (int ch = 0; ch < im.c; ++ch)
            {
                REQUIRE(res(i, j, ch) == sumWeights * (ch + 1));
            }
        }
    }

}