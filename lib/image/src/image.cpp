#include <image/image.hpp>

namespace sight
{
    template<>
    static int Image<uint8_t>::GetOpenCVType(const Image<uint8_t>& im)
    {
        return CV_8UC(im.c);
    }

} // namespace sight