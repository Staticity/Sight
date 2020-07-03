#pragma once

#include <features/brief_pattern.hpp>
#include <features/fast.hpp>
#include <features/feature.hpp>
#include <features/harris.hpp>
#include <image/image_ops.hpp>

namespace sight
{
    struct OrbDescriptor
    {
        const static int NUM_BYTES = sight::ORB_NUM_BITS / 8;

        OrbDescriptor()
            : feat()
        {
            memset(bits, 0, sizeof(bits));
        }

        OrbDescriptor(const Feature& f)
            : feat(f)
        {
            memset(bits, 0, sizeof(bits));
        }

        void setBit(int i, bool v)
        {
            assert(i >= 0 && i < ORB_NUM_BITS);

            const int j = i / 8;
            const int b = i - j * 8; // const b = i % 8;

            if (v)
            {
                bits[j] |= (1 << b);
            }
            else
            {
                bits[j] &= ~(1 << b);
            }
        }

        bool operator<(const OrbDescriptor& d) const
        {
            if (feat != d.feat)
                return feat < d.feat;
            
            for (int i = 0; i < NUM_BYTES; ++i)
            {
                if (bits[i] != d.bits[i])
                    return bits[i] < d.bits[i];
            }

            return false;
        }

        bool operator>(const OrbDescriptor& d) const
        {
            if (feat != d.feat)
                return feat > d.feat;
            
            for (int i = 0; i < NUM_BYTES; ++i)
            {
                if (bits[i] != d.bits[i])
                    return bits[i] > d.bits[i];
            }

            return false;
        }

        int distance(const OrbDescriptor& d) const
        {
            int dist = 0;

            // Compute the hamming distance between the
            // binary descriptors.
            for (int i = 0; i < NUM_BYTES; ++i)
            {
                uint8_t delta = bits[i] ^ d.bits[i];
                while (delta)
                {
                    dist += delta & 1;
                    delta >>= 1;
                }
            }

            return dist;
        }

        Feature feat;
        uint8_t bits[NUM_BYTES];
    };

    template <typename T>
    void GenerateBinaryDescriptor(
        const Image<T>& im,
        OrbDescriptor& d)
    {
        // Round the keypoint coordinates to the nearest pixel in its
        // corresponding level.
        const auto fx = d.feat.x;
        const auto fy = d.feat.y;
        const int x = int(round(fx));
        const int y = int(round(fy));

        // Pre-compute cos/sin for warping the 
        const float cosp = cos(d.feat.angle);
        const float sinp = sin(d.feat.angle);

        // Iterate over all of the pattern's offsets
        auto ptr = &ORB_BIT_PATTERN_31x31Patch[0];
        for (int i = 0; i < OrbDescriptor::NUM_BYTES; ++i)
        {
            // Directly compute each bit for the descriptor
            uint8_t bit = 0;
            for (int j = 0; j < 8; ++j, ptr += 4)
            {
                // Initial Position 1
                const int ix0 = ptr[0];
                const int iy0 = ptr[1];
                
                // Initial Position 2
                const int ix1 = ptr[2];
                const int iy1 = ptr[3];

                // Transform the coordinates by the angle of the corner.
#if 0
                // Clamp by the image's boundaries to generate a repeating
                // border for any out of bounds pixels.
                const int x0 = std::clamp(int(round(fx + cosp * ix0 - sinp * iy0)), 0, im.w - 1);
                const int y0 = std::clamp(int(round(fy + sinp * ix0 + cosp * iy0)), 0, im.h - 1);

                const int x1 = std::clamp(int(round(fx + cosp * ix1 - sinp * iy1)), 0, im.w - 1);
                const int y1 = std::clamp(int(round(fy + sinp * ix1 + cosp * iy1)), 0, im.h - 1);
#else
                // Assume border-problem points have been filtered out.
                const int x0 = int(round(fx + cosp * ix0 - sinp * iy0));
                const int y0 = int(round(fy + sinp * ix0 + cosp * iy0));

                const int x1 = int(round(fx + cosp * ix1 - sinp * iy1));
                const int y1 = int(round(fy + sinp * ix1 + cosp * iy1));
#endif

                // TODO: Determine whether we should use the original feature's
                // coordinates, or the rounded ones. Technically, the original
                // will be the most 'accurate'.

                // Look up the pixel values
                const auto v0 = im(y0, x0);
                const auto v1 = im(y1, x1);

                // Store the comparison result into the byte
                bit |= ((v0 < v1) << j);
            }

            // Store the byte into the descriptor
            d.bits[i] = bit;
        }
    }

    template <typename T>
    std::vector<OrbDescriptor> ComputeORB(
        const Pyramid<T>& pyr,
        const std::vector<Feature>& feats,
        const int nMaxFeatures = -1)
    {
        std::vector<OrbDescriptor> descs(feats.size());

        // Normalizer
        const int harrisR = 3;
        float scale = 1.f / (4 * (harrisR * 2 + 1) * 255.0f);
        float scalesq = pow(scale, 2);

        // Compute the necessary padding to perform all of the
        // image operations.
        const int briefR = 31;
        const int gaussR = 3;
        const int padding = std::max({briefR, gaussR, harrisR});

        std::vector<Image<float>> paddedLevels(pyr.NumLevels());
        for (int i = 0; i < pyr.NumLevels(); ++i)
        {
            const auto& im = pyr.GetLevel(i);

            // Compute a blurred, padded copy of the pyramid level
            const auto blur = GaussianBlur<T, float>(PadView(im, gaussR), 2.0f, gaussR);
            paddedLevels[i] = PadView(blur, padding);
        }

        // Score all of our features to select only the best ones.
        for (int i = 0; i < feats.size(); ++i)
        {
            const auto& im = paddedLevels[feats[i].level];

            // Copy the feature into the descriptor
            descs[i].feat = feats[i];

            const int cx = int(round(feats[i].x));
            const int cy = int(round(feats[i].y));

            // Compute the harris corner score for the corner's weight
            const float harrisK = 0.04f; // Harris corner constant. Usually chosen between [.04, .06]
            descs[i].feat.response = FiniteHarrisScore(im, cx, cy, harrisR, harrisK);

            // We can normalize the responses beteween 0 and 1
            descs[i].feat.response *= scalesq;
        }

        // Keep only the best N features
        std::sort(descs.begin(), descs.end(), std::greater<OrbDescriptor>());
        if (nMaxFeatures > 0)
        {
            descs.resize(std::min(descs.size(), size_t(nMaxFeatures)));
        }

        // Generate descriptors for the top features
        for (auto& desc : descs)
        {
            const auto& im = paddedLevels[desc.feat.level];

            // Compute the angle based on FAST's central moment
            desc.feat.angle = AngleFAST(im, desc.feat.x, desc.feat.y);

            GenerateBinaryDescriptor(im, desc);
        }

        return descs;
    }
}