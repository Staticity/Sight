#include <image/image.hpp>

namespace sight
{
    template<>
    static int Image<uint8_t>::GetOpenCVType(const Image<uint8_t>& im)
    {
        return CV_8UC(im.c);
    }

    template<>
    static int Image<uint16_t>::GetOpenCVType(const Image<uint16_t>& im)
    {
        return CV_16UC(im.c);
    }

    template<>
    static int Image<int>::GetOpenCVType(const Image<int>& im)
    {
        return CV_32SC(im.c);
    }

    template<>
    static int Image<float>::GetOpenCVType(const Image<float>& im)
    {
        return CV_32FC(im.c);
    }

    template<>
    static int Image<double>::GetOpenCVType(const Image<double>& im)
    {
        return CV_64FC(im.c);
    }


} // namespace sight