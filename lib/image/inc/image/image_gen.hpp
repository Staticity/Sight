#pragma once

#include <image/image.hpp>
#include <Eigen/Eigen>

#define _USE_MATH_DEFINES
#include <math.h>

namespace sight
{
    template <typename T, typename S>
    void CopyTo(const Image<T>& in, Image<S>& out)
    {
        assert(in.w == out.w);
        assert(in.h == out.h);

        for (int i = 0; i < out.h; ++i)
        {
            for (int j = 0; j < out.w; ++j)
            {
                for (int ch = 0; ch < out.c; ++ch)
                {
                    out(i, j, ch) = in(i, j, ch % in.c);
                }
            }
        }
    }

    template <typename T>
    Image<T> MergeHorizontal(const Image<T>& l, const Image<T>& r)
    {
        const int w = l.w + r.w;
        const int h = std::max(l.h, r.h);
        const int c = std::max(l.c, r.c);

        Image<T> res(w, h, c);
        CopyTo(l, res(l.GetROI()));
        CopyTo(r, res(r.GetROI().translate(l.w, 0)));
        
        return res;
    }

    template <typename T>
    Image<T> MergeVertical(const Image<T>& u, const Image<T>& d)
    {
        const int w = std::max(u.w, d.w);
        const int h = u.h + d.h;
        const int c = std::max(u.c, d.c);

        Image<T> res(w, h, c);
        CopyTo(u, res(u.GetROI()));
        CopyTo(d, res(d.GetROI().translate(u.h, 0)));

        return res;
    }

    template <typename T>
    Image<T> PadView(
        const Image<T>& im,
        const int padding)
    {
        return Pad(im, padding)(im.GetROI() + padding);
    }

    template <typename T>
    Image<T> Pad(
        const Image<T>& im,
        const int padding)
    {
        const int p = padding;
        const int p2 = p * 2;
        Image<T> padded(im.w + p2, im.h + p2, im.c);
        const int pc = p * padded.col_step;
        
        // First, fill corners by repeating the corner element
        const auto s00 = im.at(0, 0);
        const auto s01 = im.at(0, im.w - 1);
        const auto s10 = im.at(im.h - 1, 0);
        const auto s11 = im.at(im.h - 1, im.w - 1);

        for (int i = 0; i < p; ++i)
        {
            // The top left corners of all corners to fill in padding
            auto d00 = padded.at(i, 0);
            auto d01 = padded.at(i, im.w + p);
            auto d10 = padded.at(im.h + p + i, 0);
            auto d11 = padded.at(im.h + p + i, im.w + p);

            for (int j = 0; j < p; ++j)
            {
                for (int c = 0; c < im.c; ++c)
                {
                    *(d00++) = s00[c];
                    *(d01++) = s01[c];
                    *(d10++) = s10[c];
                    *(d11++) = s11[c];
                }
            }
        }

        // Now fill in all of the repeated edges

        // Fill in the top and bottom first
        auto top = im.row(0);
        auto bot = im.row(im.h - 1);

        const int nrow = im.w * im.col_step;
        for (int i = 0; i < p; ++i)
        {
            auto topdst = padded.at(i, p);
            auto botdst = padded.at(im.h + p + i, p);

            for (int j = 0; j < nrow; ++j)
            {
                topdst[j] = top[j];
                botdst[j] = bot[j];
            }
        }

        // Now we'll fill in all of the central
        // region in one go!
        for (int i = 0; i < im.h; ++i)
        {
            auto dst = padded.at(p + i, 0);
            auto src = im.row(i);

            // Repeat the left edge
            for (int j = 0; j < p; ++j)
            {
                for (int k = 0; k < im.c; ++k)
                {
                    *(dst++) = src[k];
                }
            }

            // Copy the original image
            for (int j = 0; j < nrow; ++j)
            {
                *(dst++) = src[j];
            }

            // Repeat the right edge
            for (int j = 0; j < p; ++j)
            {
                for (int k = 0; k < im.c; ++k)
                {
                    *(dst++) = src[nrow - im.c + k];
                }
            }
        }

        return padded;
    }

    template <typename T>
    Image<T> CreatePaddedImage(int w, int h, int c, int pad)
    {
        const int p2 = pad + pad;
        Image<T> im(w + p2, h + p2, c);
        return PadView(im, pad);
    }

    // Creates a grayscale corner whose corner
    // forms an interior angle of theta.
    template <typename T>
    Image<T> CreateCorner(int size, double theta, const T& b, const T& w)
    {
        Image<T> im(size, size, 1);

        // Clamp theta
        theta = std::clamp(theta, 0.0, 2 * M_PI);

        const double cx = size / 2.0;
        const double cy = cx;

        for (int i = 0; i < im.h; ++i)
        {
            for (int j = 0; j < im.w; ++j)
            {
                const double y = i - cy;
                const double x = j - cx;

                // We negate y since image coords have
                // y flipped compared to euclidean coords.
                double r = atan2(-y, x);
                if (r < 0)
                {
                    r += 2 * M_PI;
                }
                im(i, j) = (r < theta) ? w : b;
            }
        }

        return im;
    }

} // namespace sight