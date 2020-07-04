#pragma once

#include <image/image.hpp>
#include <Eigen/Core>
#include <type_traits>

#define _USE_MATH_DEFINES
#include <math.h>

namespace sight
{
    template <typename Out, typename T, typename S>
    inline Out ScalePx(const T x, const S s)
    {
        return static_cast<Out>((x + .5f) * s - .5f);
    }

    template <typename Out, typename T, typename S = float>
    inline Out BilinearInterpolate(
        const T* ptr,
        S unitx,
        S unity,
        int row_step,
        int col_step)
    {
        // Bilinear Interpolation will look at the 4 pixels bordering
        // a particular sub-pixel and linearly interpolate the 4 corners
        // based on their distance to each one.

        // 6 mults, 3 adds, 3 subs
        return static_cast<Out>(
            (1 - unitx) * (ptr[    0   ] * (1 - unity) + ptr[row_step           ] * unity) +
            unitx       * (ptr[col_step] * (1 - unity) + ptr[row_step + col_step] * unity));
    }
 
    template <typename Out, typename T, typename S = float>
    inline Out BilinearInterpolate(const Image<T>& im, S x, S y, int ch = 1)
    {
        const int x0 = int(x);
        const int y0 = int(y);

        return BilinearInterpolate<Out, T, S>(
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
    
    template <typename Out, typename T>
    Image<Out> Sobel3x3_Horizontal(const Image<T>& im)
    {
        const std::vector<int> kernelH = {-1, 0, 1};
        const std::vector<int> kernelV = { 1, 2, 1};
        return ConvolveSeparable<T, int, Out>(im, kernelH, kernelV);
    }
    
    template <typename Out, typename T>
    Image<Out> Sobel3x3_Vertical(const Image<T>& im)
    {
        const std::vector<int> kernelH = { 1, 2, 1};
        const std::vector<int> kernelV = {-1, 0, 1};
        return ConvolveSeparable<T, int, Out>(im, kernelH, kernelV);
    }

    template <typename S>
    std::vector<S> GaussianKernel(int radius, S sigma)
    {
        const int size = radius * 2 + 1;
        std::vector<S> kernel(size);

        const S twosigmasq = S(2) * sigma * sigma;
        const S scale = S(1) / (twosigmasq * S(M_PI));

        S sum = S(0);
        auto kern = kernel.begin() + radius;
        for (int i = -radius; i <= radius; ++i)
        {
            const S v = exp(-(i * i) / twosigmasq) * scale;
            kern[i] = v;
            sum += v;
        }

        const S normalizer = S(1) / sum;
        for (auto& v : kernel)
        {
            v *= normalizer;
        }

        return kernel;
    }

    template <typename T, typename Out>
    Image<Out> GaussianBlur(
        const Image<T>& im,
        const float sigma,
        const int radius)
    {
        const std::vector<float> kernel = GaussianKernel(radius, sigma);
        return ConvolveSeparable<T, float, Out>(im, kernel, kernel);
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
        const std::vector<S>& kernelH,
        const std::vector<S>& kernelV)
    {
        assert(kernelH.size() % 2 == 1);
        assert(kernelV.size() % 2 == 1);

        const int padding = int(std::max(kernelH.size() / 2, kernelV.size() / 2));
        const auto tmp = ConvolveHorizontal1D<T, S, Out>(im, kernelH);
        const auto view = Pad(tmp, padding)(Roi(0, 0, im.w, im.h) + padding);
        const auto res = ConvolveVertical1D<Out, S, Out>(view, kernelV);

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

        auto srcRow = im.begin();
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