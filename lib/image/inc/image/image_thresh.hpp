#pragma once

#include <type_traits>

#include <image/image.hpp>
#include <image/integral_image.hpp>

namespace sight
{

    template <typename T>
    Image<uint8_t> AdaptiveThreshold(
        const Image<T>& im,
        const int radius,
        const T bias,
        const uint8_t falseVal = 0,
        const uint8_t trueVal = 255)
    {
        using S = std::common_type<T, int>::type;

        Image<uint8_t> thresh(im.w, im.h, im.c);
        
        const int twoR = radius * 2;
        const int n = (twoR + 1) * (twoR + 1);

        IntegralImage<S> block(im);
        for (int i = 0; i < im.h; ++i)
        {
            auto src = im.row(i);
            auto dst = thresh.row(i);
            for (int j = 0; j < im.w; ++j)
            {
                const T s = block.sum(i - radius, j - radius, twoR, twoR);

                // s/n - c > v
                if (s > (src[j] + bias) * n)
                {
                    dst[j] = trueVal;
                }
                else
                {
                    dst[j] = falseVal;
                }
            }
        }

        return thresh;
    }
}