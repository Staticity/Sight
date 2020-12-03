#pragma once

#include <linear/vec.hpp>

#include <string>
#include <algorithm>
#include <memory>

#define _USE_MATH_DEFINES
#include <math.h>

namespace sight
{

    template <typename S>
    class CameraModel
    {
    public:

        virtual ~CameraModel() = default;

        virtual bool Project(S x, S y, S z, S& u, S& v, S* Jp = 0, S* Jxyz = 0) const = 0;

        virtual void Project(const Vec2<S>& xzyz, Vec2<S>& uv, S* Jp = 0, S* Jxyz = 0)
        {
            Project(xzyz(0), xzyz(1), S(1), uv(0), uv(1), Jp, Jxyz);
        }

        virtual void Project(const Vec3<S>& xyz, Vec2<S>& uv, S* Jp = 0, S* Jxyz = 0)
        {
            Project(xyz(0), xyz(1), xyz(2), uv(0), uv(1), Jp, Jxyz);
        }

        virtual bool Unproject(S u, S v, S& xz, S& yz) const = 0;

        virtual void Unproject(const Vec2<S>& uv, Vec2<S>& xzyz)
        {
            Unproject(uv(0), uv(1), xzyz(0), xzyz(1));
        }

        virtual void Unproject(const Vec2<S>& uv, Vec3<S>& xzyz)
        {
            Unproject(uv(0), uv(1), xzyz(0), xzyz(1));
            xzyz(2) = S(1);
        }

        virtual bool JxyzToJxzyz(S z, const S* Jxyz, S* Jxzyz) const
        {
            // Simply, we can exclude the last column of Jxyz and multiply everything else by z.
            Jxzyz[0] = z * Jxyz[0];
            Jxzyz[1] = z * Jxyz[1];
            Jxzyz[2] = z * Jxyz[0 + 3];
            Jxzyz[3] = z * Jxyz[1 + 3];
            return true;
        }

        virtual S& Param(int i) = 0;
        virtual const S& Param(int i) const = 0;
        virtual int NumParams() const = 0;

        virtual bool LoadModel(const CameraModel& model)
        {
            if (model.Name() != Name())
            {
                return false;
            }

            for (int i = 0; i < model.NumParams(); ++i)
            {
                Param(i) = model.Param(i);
            }
            return true;
        }

        virtual std::string Name() const = 0;
        virtual std::unique_ptr<CameraModel> Clone() const = 0;

        virtual S ComputeProjectiveRadius(
            const float percentile = .9,
            const S maxRadius = S(4.0)) const
        {
            // Generate a set of rays, symmetrically radiating out from
            // [0, 0, 1] to see where unprojection begins to fail. We'll
            // take the value which gives us at least `percentile` valid
            // projections.
            const int N = 90;
            const double thetaStep = (2 * M_PI) / (N - 1);
            const S maxRadiusSq = maxRadius * maxRadius;

            // Determines how much fast we'll traverse across the plane
            // at z=1. The bigger it is, the faster the search, but less
            // accurate. We could do a binary search for arbitrary precision,
            // but it's not gauranteed that the projection is valid
            // and then strictly fails beyond that critical point.
            const S metricStep = S(.1);

            S radii[N];

            for (int i = 0; i < N; ++i)
            {
                const S theta = S(thetaStep * i);
                const S dx = cos(theta) * metricStep;
                const S dy = sin(theta) * metricStep;

                S x = dx;
                S y = dy;
                const S z = S(1);

                S u, v;
                S xz, yz; // unused
                while (Project(x, y, z, u, v) && Unproject(u, v, xz, yz))
                {
                    x += dx;
                    y += dy;

                    // We've gone far enough
                    if (x * x + y * y > maxRadiusSq)
                    {
                        break;
                    }
                }

                // We remove one more step, since this is the failing step.
                x -= dx;
                y -= dy;
                radii[i] = x * x + y * y;
            }

            std::sort(radii, radii + N);
            const int idx = int(round((N - 1) * (1.0f - percentile)));
            return sqrt(radii[idx]);
        }

        virtual bool IterativeUnproject(S u, S v, S& xz, S& yz, const int maxIterations) const
        {
            // Assume xz and yz are properly initialized

            S H[4];
            S g[2];
            S Jxyz[6];
            S Jxzyz[4];

            // Run up to `maxIterations` number of gauss newton optimization steps
            for (int i = 0; i < maxIterations; ++i)
            {
                S up, vp;
                Project(xz, yz, S(1), up, vp, nullptr, Jxyz);
                JxyzToJxzyz(S(1), Jxyz, Jxzyz); // Should do nothing, since z=1!

                // Compute the residual
                const S du = up - u;
                const S dv = vp - v;

                if (du * du + dv * dv < std::numeric_limits<S>::epsilon())
                {
                    return true;
                }

                // Compute the hessian approximation J^T * J
                H[0] = Jxzyz[0] * Jxzyz[0] + Jxzyz[2] * Jxzyz[2];
                H[1] = Jxzyz[0] * Jxzyz[1] + Jxzyz[2] * Jxzyz[3];
                H[2] = H[1];
                H[3] = Jxzyz[1] * Jxzyz[1] + Jxzyz[3] * Jxzyz[3];

                // Compute J^T * r
                g[0] = Jxzyz[0] * du + Jxzyz[2] * dv;
                g[1] = Jxzyz[1] * du + Jxzyz[3] * dv;

                // Now, we need to solve Hv = -g
                const S det = H[0] * H[3] - H[1] * H[2];

                // Is H invertible?
                if (det == S(0))
                {
                    return false;
                }

                // The inverse of a 2x2 matrix is:
                //
                //  [a b]^-1  =  1/det [d -b]
                //  [c d]              [-c a]
                //
                // So, then our update will be:
                //
                // update = H^-1 * -g
                const S invdet = S(1) / det;
                xz += -invdet * (H[3] * g[0] - H[1] * g[1]);
                yz += -invdet * (H[0] * g[1] - H[2] * g[0]);
            }

            return false;
        }

    };
}