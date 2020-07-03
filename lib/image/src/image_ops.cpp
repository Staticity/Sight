#include <image/image_ops.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

namespace sight
{
    std::vector<float> GaussianKernel(int radius, float sigma)
    {
        const int size = radius * 2 + 1;
        std::vector<float> kernel(size);

        const float twosigmasq = 2 * sigma * sigma;
        const float scale = 1.0f / (twosigmasq * float(M_PI));

        float sum = 0.f;
        auto kern = kernel.begin() + radius;
        for (int i = -radius; i <= radius; ++i)
        {
            const float v = exp(-(i * i) / twosigmasq) * scale;
            kern[i] = v;
            sum += v;
        }

        const float normalizer = 1.f / sum;
        for (auto& v : kernel)
        {
            v *= normalizer;
        }

        return kernel;
    }

    Image<uint8_t> RadialMask(int r)
    {
        const int size = r * 2 + 1;
        Image<uint8_t> mask(size, size, 1);

        const int r2 = r * r;
        auto ptr = mask.begin() + r;
        for (int i = -r; i <= r; ++i, ptr += size)
        {
            const int i2 = i * i;
            for (int j = -r; j <= r; ++j)
            {
                ptr[j] = (i2 + j * j) <= r2;
            }
        }

        return mask;
    }

}