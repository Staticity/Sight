#include <image/image_ops.hpp>

namespace sight
{
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