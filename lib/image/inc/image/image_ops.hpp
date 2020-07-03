#pragma once

#include <image/image.hpp>
#include <Eigen/Core>
#include <type_traits>

namespace sight
{
    template <typename Out, typename T, typename S>
    inline Out ScalePx(const T x, const S s)
    {
        return static_cast<Out>((x + .5f) * s - .5f);
    }

    template <typename Out, typename T>
    inline Out BilinearInterpolate(
        const T* ptr,
        float unitx,
        float unity,
        int row_step,
        int col_step)
    {
        // Bilinear Interpolation will look at the 4 pixels bordering
        // a particular sub-pixel and linearly interpolate the 4 corners
        // based on their distance to each one.

        const float p00 = (1 - unitx) * (1 - unity);
        const float p01 = (1 - unitx) * unity;
        const float p10 = unitx * (1 - unity);
        const float p11 = unitx * unity;
        return static_cast<Out>(
            p00 * ptr[0] +
            p01 * ptr[row_step] +
            p10 * ptr[col_step] +
            p11 * ptr[row_step + col_step]
        );
    }
 
    template <typename Out, typename T>
    inline Out BilinearInterpolate(const Image<T>& im, float x, float y, int ch = 1)
    {
        const int x0 = int(x);
        const int y0 = int(y);

        assert(x0 >= 0);
        assert(y0 >= 0);
        assert(x0 + 1 < im.w);
        assert(y0 + 1 < im.h);

        return BilinearInterpolate<T, Out>(
            im.at(y0, x0, ch),
            x - x0,
            y - y0,
            im.row_step,
            im.col_step);
    }

    template <typename Out, typename T>
    void Resize(const Image<T>& im, Image<Out>& res)
    {
        assert(im.c == res.c);
        const float sx = im.w / float(res.w);
        const float sy = im.h / float(res.h);

        for (int i = 0; i < res.h; ++i)
        {
            const float y = ScalePx<float>(i, sy);
            
            auto dst = res.row(i);
            for (int j = 0; j < res.w; ++j)
            {
                const float x = ScalePx<float>(j, sx);

                for (int ch = 0; ch < res.c; ++ch)
                {
                    *(dst++) = BilinearInterpolate<Out>(im, x, y, ch);
                }
            }
        }
    }

    template <typename Out, typename T>
    Image<Out> Resize(const Image<T>& im, const float sx, const float sy)
    {
        const int w = int(im.w * sx);
        const int h = int(im.h * sy);
        auto res = CreatePaddedImage<Out>(w, h, im.c, 0);
        Resize(im, res);
        return res;
    }

    std::vector<float> GaussianKernel(int radius, float sigma);

    template <typename T, typename Out>
    Image<Out> GaussianBlur(
        const Image<T>& im,
        const float sigma,
        const int radius)
    {
        const std::vector<float> kernel = GaussianKernel(radius, sigma);
        return ConvolveSeparable<T, float, Out>(im, kernel);
    }

    template <typename T, typename S, typename Out>
    void ConvolveHorizontal1D(
        const Image<T>& im,
        const std::vector<S>& kernel,
        Image<Out>& res)
    {
        assert(kernel.size() % 2 == 1);

        // Center the kernel pointer
        const int r = int(kernel.size() / 2);
        const auto kern = kernel.begin() + r;

        // Convolve horizontally
        using TS = std::common_type<T, S>::type;
        auto dst = res.begin();
        for (int i = 0; i < res.h; ++i)
        {
            const auto row = im.row(i);
            for (int j = 0; j < res.w; ++j)
            {
                const auto px = row + j * im.col_step;
                for (int ch = 0; ch < res.c; ++ch)
                {
                    const auto src = px + ch;
                    
                    // Convolve src into dst
                    TS v = TS(0);
                    for (int k = -r; k <= r; ++k)
                    {
                        v += src[k * im.col_step] * kern[k];
                    }
                    *(dst++) = static_cast<Out>(v);

                }
            }
        }
    }
    
    template <typename T, typename S, typename Out>
    Image<Out> ConvolveHorizontal1D(
        const Image<T>& im,
        const std::vector<S>& kernel)
    {
        Image<Out> res(im.w, im.h, im.c);
        ConvolveHorizontal1D(im, kernel, res);
        return res;
    }
    
    template <typename T, typename S, typename Out>
    void ConvolveVertical1D(
        const Image<T>& im,
        const std::vector<S>& kernel,
        Image<Out>& res)
    {
        assert(kernel.size() % 2 == 1);

        // Center the kernel pointer
        const int r = int(kernel.size() / 2);
        const auto kern = kernel.begin() + r;

        // Convolve horizontally
        using TS = std::common_type<T, S>::type;
        auto dst = res.begin();
        for (int i = 0; i < res.h; ++i)
        {
            const auto row = im.row(i);
            for (int j = 0; j < res.w; ++j)
            {
                const auto px = row + j * im.col_step;
                for (int ch = 0; ch < res.c; ++ch)
                {
                    const auto src = px + ch;
                    
                    // Convolve src into dst
                    TS v = TS(0);
                    for (int k = -r; k <= r; ++k)
                    {
                        v += src[k * im.row_step] * kern[k];
                    }
                    *(dst++) = static_cast<Out>(v);

                }
            }
        }
    }

    template <typename T, typename S, typename Out>
    Image<Out> ConvolveVertical1D(
        const Image<T>& im,
        const std::vector<S>& kernel)
    {
        Image<Out> res(im.w, im.h, im.c);
        ConvolveVertical1D(im, kernel, res);
        return res;
    }

    template <typename T, typename S, typename Out>
    Image<Out> ConvolveSeparable(
        const Image<T>& im,
        const std::vector<S>& kernel)
    {
        assert(kernel.size() % 2 == 1);

        const int padding = int(kernel.size() / 2);
        const auto tmp = ConvolveHorizontal1D<T, S, Out>(im, kernel);
        const auto view = Pad(tmp, padding)(Roi(0, 0, im.w, im.h) + padding);
        const auto res = ConvolveVertical1D<Out, S, Out>(view, kernel);

        return res;
    }

    Image<uint8_t> RadialMask(int r);

    template <typename T>
    Eigen::Vector2f CentralMomentWithMask(
        const Image<T>& im,
        int x,
        int y,
        const Image<uint8_t>& mask)
    {
        assert(im.c == 1);
        assert(mask.c == 1);

        // Assumption: The mask is square and of odd size
        assert(mask.w == mask.h);
        assert(mask.w % 2 == 1);

        float m00 = 0.0;
        float m01 = 0.0;
        float m10 = 0.0;
        
        const int r = mask.w / 2;
        auto ptr = im.row(y - r) + x;
        auto maskptr = mask.begin() + r;
        for (int i = -r; i <= r; ++i)
        {
            for (int j = -r; j <= r; ++j)
            {
                if (maskptr[j])
                {
                    const auto& v = ptr[j];
                    m00 += v;
                    m01 += i * v;
                    m10 += j * v;
                }
            }

            ptr += im.row_step;
            maskptr += mask.row_step;
        }

        return Eigen::Vector2f(x + m10 / m00, y + m01 / m00);
    }


    template <typename T>
    Eigen::Vector2f CentralMoment(
        const Image<T>& im,
        int x,
        int y,
        int r)
    {
        assert(im.c == 1);
        const int r2 = r * r;

        float m00 = 0.0;
        float m01 = 0.0;
        float m10 = 0.0;

        auto ptr = im.row(y - r) + x;
        for (int i = -r; i <= r; ++i, ptr += im.row_step)
        {
            const int i2 = i * i;
            for (int j = -r; j <= r; ++j)
            {
                if (i2 + j * j <= r2)
                {
                    const auto& v = ptr[j];
                    m00 += v;
                    m01 += i * v;
                    m10 += j * v;
                }
            }
        }

        return Eigen::Vector2f(x + m10 / m00, y + m01 / m00);
    }

    template <typename T>
    void Fill(Image<T>& im, const T& v)
    {
        for (auto it = im.begin(); it < im.end(); ++it)
        {
            *it = v;
        }
    }

    template <typename T>
    Image<T> ToGrayscale(const Image<T>& im)
    {
        if (im.c == 1) return im;
        Image<T> gray(im.w, im.h, 1);

        assert(im.c == 3);
        if (im.c != 3) return gray;

        // Linear coefficients for determining
        // perceieved intensity.
        const float r = 0.2989f;
        const float g = 0.5870f;
        const float b = 0.1140f;

        // Store grayscale results into `gray`.
        const T* srcRow = im.begin();
        T* dstRow = gray.begin();
        int n = im.w * im.h;
        for (int i = 0; i < im.h; ++i)
        {
            for (int j = 0, k = 0; j < im.w; ++j)
            {
                dstRow[j] = T(srcRow[k++] * b + srcRow[k++] * g + srcRow[k++] * r);
            }
            srcRow += im.row_step;
            dstRow += gray.row_step;
        }

        return gray;
    }

    template <typename T>
    Image<T> ToRGB(const Image<T>& im)
    {
        if (im.c == 3) return im;
        assert(im.c == 1);

        Image<T> rgb(im.w, im.h, 3);

        const auto srcRow = im.begin();
        auto dstRow = rgb.begin();
        for (int i = 0; i < im.h; ++i)
        {
            for (int j = 0, k = 0; j < im.w; ++j)
            {
                const auto v = srcRow[j];
                dstRow[k++] = v;
                dstRow[k++] = v;
                dstRow[k++] = v;
            }
            srcRow += im.row_step;
            dstRow += rgb.row_step;
        }

        return rgb;
    }
    
} // namespace sight