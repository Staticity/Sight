#pragma once

#include <array>
#include <calibration/cameramodel.hpp>

namespace sight
{
    template <typename S>
    class PinholeModel : public ICameraModel<S>
    {
    public:
        enum
        {
            FX = 0,
            FY,
            CX,
            CY,
            SKEW,
            NUM_PARAMS
        };

        PinholeModel(
            const S fx = S(1),
            const S fy = S(1),
            const S cx = S(0),
            const S cy = S(0),
            const S skew = S(0))
        {
            p[0] = fx;
            p[1] = fy;
            p[2] = cx;
            p[3] = cy;
            p[4] = skew;
        }

        bool Project(
            S x, S y, S z,
            S& u, S& v,
            S* Jp = nullptr,
            S* Jxyz = nullptr) const override
        {
            // [fx sk cx] [x]   [fx*x + sk*y + cx]
            // [ 0 fy cy] [y] = [fy*y + cy]
            // [ 0  0  1] [z]   [z]
            //
            // Then perform perspective division by dividing
            // by z.
            const S invz = S(1) / z;
            const S xp = x * invz;
            const S yp = y * invz;

            u = p[FX] * xp + p[SKEW] * yp + p[CX];
            v = p[FY] * yp + p[CY];

            if (Jp)
            {
                S* Ju = Jp;
                S* Jv = Jp + NUM_PARAMS;
    
                Ju[FX] = xp;
                Jv[FX] = S(0);

                Ju[FY] = S(0);
                Jv[FY] = yp;

                Ju[CX] = S(1);
                Jv[CX] = S(0);

                Ju[CY] = S(0);
                Jv[CY] = S(1);

                Ju[SKEW] = yp;
                Jv[SKEW] = S(0);
            }

            if (Jxyz)
            {
                S* Ju = Jxyz;
                S* Jv = Jxyz + 3;

                Ju[0] = p[FX] * invz;
                Jv[0] = S(0);

                Ju[1] = p[SKEW] * invz;
                Jv[1] = p[FY] * invz;

                Ju[2] = -(p[FX] * xp + p[SKEW] * yp) * invz;
                Jv[2] = -(p[FY] * yp) * invz;
            }

            return true;
        }

        bool Unproject(S u, S v, S& xz, S& yz) const override
        {
            // Need to compute:
            //
            // [fx sk cx]^-1  [u]
            // [ 0 fy cy]     [v]
            // [ 0  0  1]     [1]
            //
            // Unprojecting the point onto the plane z=1.

            // We have:
            //
            //     v = fy*(y/z) + cy
            //
            // So then with z=1, we have:
            //
            yz = (v - p[CY]) / p[FY];

            // With y solved, we now have:
            //
            //     u = (fx*x + sk*y)/z + cx
            //
            // So then with z=1, we have:
            //
            xz = (u - p[CX] - p[SKEW] * yz) / p[FX];

            return true;
        }

        inline S& Param(int i) override { return p[i]; }
        inline const S& Param(int i) const override  { return p[i]; }
        inline int NumParams() const override { return NUM_PARAMS; }

        static std::string ModelName() { return "pinhole"; }
        inline std::string Name() const override { return ModelName(); }

        std::unique_ptr<ICameraModel<S>> Clone() const override
        {
            std::unique_ptr<ICameraModel<S>> clone(new PinholeModel<S>());
            clone->LoadModel(*this);
            return clone;
        }

        std::array<S, NUM_PARAMS> p;
    };
}
