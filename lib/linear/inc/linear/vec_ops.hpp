#pragma once

#include <linear/vec.hpp>

namespace sight
{
    template <int N, typename S>
    S Dot(const S* u, const S* v)
    {
        S d = S(0);
        for (int i = 0; i < N; ++i)
        {
            d += u[i] * v[i];
        }
        return d;
    }

    template <typename S>
    void Cross3(const S* u, const S* v, S* uxv)
    {
        // [  0 -uz  uy] [vx]   [uy*vz - uz*vy]
        // [ uz   0 -ux] [vy] = [uz*vx - ux*vz]
        // [-uy  ux   0] [vz]   [ux*vy - uy*vx]
        constexpr int x = 0;
        constexpr int y = 1;
        constexpr int z = 2;
        uxv[x] = (u[y] * v[z]) - (u[z] * v[y]);
        uxv[y] = (u[z] * v[x]) - (u[x] * v[z]);
        uxv[z] = (u[x] * v[y]) - (u[y] * v[x]);
    }
}
