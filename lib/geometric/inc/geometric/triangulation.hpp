#pragma once

#include <linear/vec.hpp>
#include <lie/se3.hpp>

namespace sight
{

    template <typename S>
    bool Compute3DPoint(
        const Vec3<S>& m0,
        const Vec3<S>& m1,
        const Vec3<S>& t,
        Vec3<S>& point)
    {
        const Vec3<S> z = m1.Cross(m0);
        const S z2 = z.Dot(z);
        const S y0 = (z.Dot(t.Cross(m1))) / z2;
        const S y1 = (z.Dot(t.Cross(m0))) / z2;

        point = m0 * y0;

        return y0 > S(0) && y1 > S(0);
    }

    // Beautiful triangulation method - simple & fast:
    // https://arxiv.org/pdf/1903.09115.pdf
    template <typename S>
    bool TriangulateL1Angular(
        const Vec3<S>& ray0,
        const Vec3<S>& ray1,
        const SE3<S>& Cam0FromCam1,
        Vec3<S>& pointInCam0)
    {
        // The paper proposes finding the optimal adjustments
        // to the unprojected rays to minimize the L1 error.
        //
        // Surprisingly, we can minimize the L1 angular error
        // by only modifying one of the rays, but not both.
        //
        // Thus, we should determine whether it's closer to project
        // a point from ray0 onto the plane formed by the
        // translation between the cameras and the ray1, or vice versa.

        const auto& x0 = ray0;
        const auto Rx1 = Cam0FromCam1.R * ray1;

        const auto& t = Cam0FromCam1.t;

        Vec3<S> m0p, m1p;
        
        // Is it cheaper to project Rx1 or x0?
        if (x0.Normalized().Cross(t).SquaredNorm() <= Rx1.Normalized().Cross(t).SquaredNorm())
        {
            // Project Rx0
            const Vec3<S> n1 = Rx1.Cross(t).Normalized();
            m0p = x0 - n1 * (x0.Dot(n1));
            m1p = Rx1;
        }
        else
        {
            // Project x1
            const Vec3<S> n0 = x0.Cross(t).Normalized();            
            m0p = x0;
            m1p = Rx1 - n0 * (Rx1.Dot(n0));
        }

        return Compute3DPoint(m0p, m1p, t, pointInCam0);
    }
}