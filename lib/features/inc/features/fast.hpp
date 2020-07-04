#pragma once

#include <image/image.hpp>
#include <image/pyramid.hpp>
#include <image/image_ops.hpp>
#include <array>
#include <vector>
#include <features/feature.hpp>

namespace sight
{
    // This is hard-coded, but we can parameterize it later.
    const int fast_radius = 3;
    const int fast_pattern = 16;
    const int seg_size = 9;
    const auto fast_mask = sight::RadialMask(fast_radius);

    const int fast_thresh = 20;

    // The score of a FAST corner is defined as
    // the largest threshold we can use to keep
    // the corner valid.
    int ScoreFAST(
        const uint8_t* corner,
        int thresh,
        const std::vector<int>& offsets)
    {
        assert(offsets.size() == 16 + 9);

        const size_t N = offsets.size();
        const auto p = *corner;

        // Pre-compute all offsets to the corner
        std::vector<int> d(N);
        for (int i = 0; i < N; ++i)
        {
            d[i] = p - corner[offsets[i]];
        }

        // How many elements do we need to check?
        const size_t K = N - seg_size;

        // Check upper bound first
        int lo = thresh;
        for (int k = 0; k < K; ++k)
        {
            int x = std::min(d[k], d[k + 1]);
            if (x <= lo) continue;
            x = std::min(x, d[k + 2]);
            x = std::min(x, d[k + 3]);
            x = std::min(x, d[k + 4]);
            x = std::min(x, d[k + 5]);
            x = std::min(x, d[k + 6]);
            x = std::min(x, d[k + 7]);
            x = std::min(x, d[k + 8]);
            lo = std::max(lo, x);
        }

        // Now the lower bound with the previous upper
        // bound as a limit.
        int hi = -lo;
        for (int k = 0; k < K; ++k)
        {
            int x = std::max(d[k], d[k + 1]);
            if (x >= hi) continue;
            x = std::max(x, d[k + 2]);
            x = std::max(x, d[k + 3]);
            x = std::max(x, d[k + 4]);
            x = std::max(x, d[k + 5]);
            x = std::max(x, d[k + 6]);
            x = std::max(x, d[k + 7]);
            x = std::max(x, d[k + 8]);
            hi = std::min(hi, x);
        }

        return -hi - 1;
    }
    
    template <typename T>
    float AngleFAST(
        const Image<T>& im,
        const float px,
        const float py)
    {
        const int x = int(round(px));
        const int y = int(round(py));
        const auto center = CentralMomentWithMask(im, x, y, fast_mask);
        return atan2(center(1) - py, center(0) - px);
    }
    
    void FindFAST(
        const Image<uint8_t>& im,
        std::vector<Feature>& feats,
        const bool useNms = true)
    {
        if (im.c != 1) return;

        // create pointer offsets for easy circular
        // neighbor accessing.
        const int N = 16 + seg_size;
        std::vector<int> pixel =
        {
             0 + -3 * im.h, // ( 0,  3)
             1 + -3 * im.h, // ( 1,  3)
             2 + -2 * im.h, // ( 2,  2)
             3 + -1 * im.h, // ( 3,  1)
             3 + 0,         // ( 3,  0)
             3 + 1 * im.h,  // ( 3, -1)
             2 + 2 * im.h,  // ( 2, -2)
             1 + 3 * im.h,  // ( 1, -3)
             0 + 3 * im.h,  // ( 0, -3)
            -1 + 3 * im.h,  // (-1, -3)
            -2 + 2 * im.h,  // (-2, -2)
            -3 + 1 * im.h,  // (-3, -1)
            -3 + 0,         // (-3,  0)
            -3 + -1 * im.h, // (-3,  1)
            -2 + -2 * im.h, // (-2,  2)
            -1 + -3 * im.h, // (-1,  3)
        };
        pixel.resize(N);

        // Add a circular portion to simplify the
        // segment test.
        for (int i = 16; i < N; ++i)
        {
            pixel[i] = pixel[i - 16];
        }

        // Classify each pixel as belonging to one of 3
        // classes:
        //
        // 0: within [p - t, p + t]
        // 1: within (-inf, p - t)
        // 2: within (p + t, inf)
        uint8_t thresh_flag[512];

        // Assume p=0, and pre-fill the bits array
        // for efficient lookups.
        for (int i = -255; i <= 255; ++i)
        {
            thresh_flag[i + 255] = (i < -fast_thresh ? 1 : (i > fast_thresh ? 2 : 0));
        }

        // Let's create buffers to be used during
        // non-maximum suppression. It will store the
        // scores and corners of all the past 3 rows.
        int corner_count[3] = {};
        std::vector<int> corner_col[3];
        std::vector<uint8_t> row_scores[3];
        for (int i = 0; i < 3; ++i)
        {
            corner_col[i].resize(im.w);
            row_scores[i].resize(im.w);
        }

        // Convenience names for the 3x3 buffers.
        auto& pprev_count = corner_count[0];
        auto& prev_count = corner_count[1];
        auto& curr_count = corner_count[2];

        auto& pprev_corners = corner_col[0];
        auto& prev_corners = corner_col[1];
        auto& curr_corners = corner_col[2];

        auto& pprev_score = row_scores[0];
        auto& prev_score = row_scores[1];
        auto& curr_score = row_scores[2];
        
        // Now iterate over every pixel and determine if it's
        // a corner. To be efficient, we'll have to manually
        // update the row/col indices.
        const uint8_t* ptr = im.row(3);

        // For every valid row..
        for (int row = 3; row < im.h - 3; ++row, ptr += 3)
        {
            // Skip the beginning of this row
            ptr += 3;

            curr_count = 0;
            memset(curr_score.data(), 0, sizeof(uint8_t) * im.w);

            // For every valid col..
            for (int col = 3; col < im.w - 3; ++col, ++ptr)
            {
                // Determine if this pixel is a corner via
                // a 9-segment test. We'll first efficiently
                // reject non-corners.
                auto p = *ptr;
                auto* flag = thresh_flag + (255 - p);

                // If pixels 0 and 8 are not above/below the threshold,
                // then there is no 9-segment.
                int f = flag[ptr[pixel[0]]] | flag[ptr[pixel[8]]];
                
                if (f == 0) continue;

                // Check the other pairs of opposites to see if
                // they also contribute the same flags.
                f &= flag[ptr[pixel[2]]] | flag[ptr[pixel[10]]];
                f &= flag[ptr[pixel[4]]] | flag[ptr[pixel[12]]];
                f &= flag[ptr[pixel[6]]] | flag[ptr[pixel[14]]];

                if (f == 0) continue;

                // Finally, include all of the pixel comparisons
                f &= flag[ptr[pixel[1]]] | flag[ptr[pixel[9]]];
                f &= flag[ptr[pixel[3]]] | flag[ptr[pixel[11]]];
                f &= flag[ptr[pixel[5]]] | flag[ptr[pixel[13]]];
                f &= flag[ptr[pixel[7]]] | flag[ptr[pixel[15]]];

                // At this point, if `f` is nonzero, then all pairs
                // of opposites share at least one of the same
                // classifications.

                // Now perform the brute-force check for corners,
                // which determines if there are 9+ neighbor
                // pixels with the same classification in a row.
                if (f & 1)
                {
                    const int th = p - fast_thresh;
                    int inarow = 0;

                    for (int k = 0; k < N; ++k)
                    {
                        if (ptr[pixel[k]] < th)
                        {
                            ++inarow;
                            if (inarow >= seg_size)
                            {
                                // Found a corner!
                                curr_corners[curr_count++] = col;
                                curr_score[col] = ScoreFAST(ptr, fast_thresh, pixel);
                                break;
                            }
                        }
                        else
                        {
                            inarow = 0;
                        }
                    }
                }

                if (f & 2)
                {
                    const int th = p + fast_thresh;
                    int inarow = 0;

                    for (int k = 0; k < N; ++k)
                    {
                        if (ptr[pixel[k]] > th)
                        {
                            ++inarow;
                            if (inarow >= seg_size)
                            {
                                // Found a corner!
                                curr_corners[curr_count++] = col;
                                curr_score[col] = ScoreFAST(ptr, fast_thresh, pixel);
                                break;
                            }
                        }
                        else
                        {
                            inarow = 0;
                        }
                    }
                }
            }

            // Now, let's add our keypoints while accounting
            // for non-maximum suppression.
            //
            // Skip the first row, because it has no previous
            // row.
            if (row != 3)
            {
                for (int k = 0; k < prev_count; ++k)
                {
                    // Consider the previous row, since it now has
                    // a `prev` and `next` buffer to compare.
                    int i = row - 1;
                    int j = prev_corners[k];
                    int s = prev_score[j];
                    if (!useNms || (
                        s > pprev_score[j - 1] && s > pprev_score[j] && s > pprev_score[j + 1] &&
                        s > prev_score[j - 1] && s > prev_score[j + 1] &&
                        s > curr_score[j - 1] && s > curr_score[j] && s > curr_score[j + 1]))
                    {
                        feats.emplace_back(float(j), float(i));
                    }
                }
            }

            // Swap the buffers
            std::swap(curr_count, pprev_count);
            std::swap(pprev_count, prev_count);

            std::swap(curr_score, pprev_score);
            std::swap(pprev_score, prev_score);

            std::swap(curr_corners, pprev_corners);
            std::swap(pprev_corners, prev_corners);
        }
    }

    std::vector<Feature> FindFAST(const Pyramid<uint8_t>& p, bool useNms = true)
    {
        std::vector<Feature> feats;
        for (int i = 0; i < p.NumLevels(); ++i)
        {
            // Store the previous size of the list of features
            size_t j = feats.size();
            FindFAST(p.GetLevel(i), feats, useNms);

            // Assign appropriate level/scale information
            // for all the newly created features
            const float sx = p.ScaleX(i);
            const float sy = p.ScaleY(i);
            for (; j < feats.size(); ++j)
            {
                feats[j].level = i;
                feats[j].scaleX = sx;
                feats[j].scaleY = sy;
            }
        }
        return feats;
    }

} // namespace sight