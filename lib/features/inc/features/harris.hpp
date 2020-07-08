#pragma once

#include <estimation/quadratic.hpp>

namespace sight
{
    enum ScoreHeuristic
    {
        HARRIS,
        SHI_TOMASI
    };

    int FiniteCornerPadding(int k)
    {
        // Normally, you only require
        // a block size of 2k + 1, but
        // we also compute finite differences
        // weighted w/ Sobel which requires
        // another padding of 1.
        return 2 * (k + 1) + 1;
    }
    
    template <typename T>
    float FiniteCornerScore(
        const Image<T>& im,
        int x,
        int y,
        int r,
        // float constantK,
        ScoreHeuristic scoreType = SHI_TOMASI)
    {
        assert(im.c == 1);
        // assert(constantK >= .04f && constantK <= .06f);

        float Ixx = 0.f;
        float Iyy = 0.f;
        float Ixy = 0.f;

        for (int i = -r; i <= r; ++i)
        {
            auto pRow = im.row(y + i) + x /* * col_step*/;
            for (int j = -r; j <= r; ++j)
            {
                // center pointer
                const auto c = pRow + j;

                // up, down, left, right
                const auto u = c[-im.row_step];
                const auto d = c[+im.row_step];
                const auto l = c[j - 1];
                const auto r = c[j + 1];

                // up left, up right, down left, down right
                const auto ul = c[-im.row_step - 1];
                const auto ur = c[-im.row_step + 1];
                const auto dl = c[+im.row_step - 1];
                const auto dr = c[+im.row_step + 1];

                // Compute I_dx and I_dy with finite differences
                // weighted by Sobel.
                const auto dx = (ul - ur) + 2 * (l - r) + (dl - dr);
                const auto dy = (ul - dl) + 2 * (u - d) + (ur - dr);

                Ixx += dx * dx;
                Iyy += dy * dy;
                Ixy += dx * dy;
            }
        }

        // Response = e1*e2 - k * (e1 + e2)^2
        //
        // where e1 and e2 are the eigen vectors of the structure
        // tensor.
        //
        // e1 * e2 is simply the determinant of the matrix, while
        // we can say the same about e1 + e2 as the trace.

        const float det = Ixx * Iyy - Ixy * Ixy;
        const float trace = Ixx + Iyy;

        if (scoreType == HARRIS)
        {
            const float constantK = .04f;
            return det - constantK * (trace * trace);
        }
        else if (scoreType == SHI_TOMASI)
        {
            Quadratic<float> q(1.f, -trace, det);
            float roots[2];
            const int nRoots = q.RealRoots(roots);
            return nRoots == 0 ? 0.f : roots[0];
        }
        else
        {
            return 0.f;
        }
    }
}