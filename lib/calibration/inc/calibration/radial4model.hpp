#pragma once

#include <cmath>
#include <calibration/cameramodel.hpp>

namespace sight
{
    template <typename S>
    class Radial4Model : public ICameraModel<S>
    {
    public:

        enum
        {
            FX = 0,
            FY,
            CX,
            CY,
            K1,
            K2,
            K3,
            K4,
            NUM_PARAMS
        };

        Radial4Model(
            S fx = S(1),
            S fy = S(1),
            S cx = S(0),
            S cy = S(0),
            S k1 = S(0),
            S k2 = S(0),
            S k3 = S(0),
            S k4 = S(0))
        {
            p[FX] = fx;
            p[FY] = fy;
            p[CX] = cx;
            p[CY] = cy;
            p[K1] = k1;
            p[K2] = k2;
            p[K3] = k3;
            p[K4] = k4;
        }

        bool Project(
            S x, S y, S z,
            S& u, S& v,
            S* Jp = nullptr,
            S* Jxyz = nullptr) const override
        {
            const S invz = S(1) / z;
            const S xp = x * invz;
            const S yp = y * invz;

            const S r2 = xp * xp + yp * yp;
            const S r1 = sqrt(r2);
            const S r3 = r2 * r1;
            const S r4 = r2 * r2;

            const S s = S(1) + p[K1] * r1 + p[K2] * r2 + p[K3] * r3 + p[K4] * r4;
            const S xpfx = xp * p[FX];
            const S ypfy = yp * p[FY];
            u = s * xpfx + p[CX];
            v = s * ypfy + p[CY];

            if (Jp)
            {
                S* Ju = Jp;
                S* Jv = Jp + NUM_PARAMS;

                Ju[FX] = s * xp;
                Jv[FX] = S(0);

                Ju[FY] = S(0);
                Jv[FY] = s * yp;

                Ju[CX] = S(1);
                Jv[CX] = S(0);

                Ju[CY] = S(0);
                Jv[CY] = S(1);

                Ju[K1] = r1 * xpfx;
                Jv[K1] = r1 * ypfy;

                Ju[K2] = r2 * xpfx;
                Jv[K2] = r2 * ypfy;

                Ju[K3] = r3 * xpfx;
                Jv[K3] = r3 * ypfy;

                Ju[K4] = r4 * xpfx;
                Jv[K4] = r4 * ypfy;
            }

            if (Jxyz)
            {
                S* Ju = Jxyz;
                S* Jv = Jxyz + 3;

                // Question: What should this be set to when the radius is 0?
                const S invr = (r1 < std::numeric_limits<S>::epsilon()) ? S(1) : S(1) / r1;
                const S ds_dr = p[K1] + (S(2) * r1 * p[K2]) + (S(3) * r2 * p[K3]) + (S(4) * r3 * p[K4]);
                const S dr_dx = xp * invz * invr;
                const S dr_dy = yp * invz * invr;
                const S dr_dz = -r1 * invz;

                const S ds_dx = ds_dr * dr_dx;
                const S ds_dy = ds_dr * dr_dy;
                const S ds_dz = ds_dr * dr_dz;

                Ju[0] = (s * p[FX] * invz) + xpfx * ds_dx;
                Jv[0] = ypfy * ds_dx;

                Ju[1] = xpfx * ds_dy;
                Jv[1] = (s * p[FY] * invz) + ypfy * ds_dy;

                const S invz2 = invz * invz;
                Ju[2] = (p[FX] * s * -x * invz2) + (p[FX] * ds_dz * xp);
                Jv[2] = (p[FY] * s * -y * invz2) + (p[FY] * ds_dz * yp);
            }

            return true;
        }

        bool Unproject(S u, S v, S& xz, S& yz) const override
        {
             // Initialize estimate with linear calibration
            xz = (u - p[CX]) / (p[FX]);
            yz = (v - p[CY]) / (p[FY]);

            return IterativeUnproject(u, v, xz, yz, 10);
        }

        S& Param(int i) override { return p[i]; };
        const S& Param(int i) const override { return p[i]; }
        int NumParams() const override { return NUM_PARAMS; }

        static std::string ModelName() { return "radial4"; }
        inline std::string Name() const override { return ModelName(); }

        std::unique_ptr<ICameraModel> Clone() const override
        {
            std::unique_ptr<ICameraModel<S>> clone(new Radial4Model<S>());
            clone->LoadModel(*this);
            return clone;
        }
        
        std::array<S, NUM_PARAMS> p;
    };
}