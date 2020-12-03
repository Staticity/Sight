#pragma once

#include <array>
#include <cmath>
#include <calibration/cameramodel.hpp>

namespace sight
{
    template <typename S>
    class EquidistantModel : public CameraModel<S>
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

        EquidistantModel(
            const S fx = S(1),
            const S fy = S(1),
            const S cx = S(0),
            const S cy = S(0),
            const S k1 = S(0),
            const S k2 = S(0),
            const S k3 = S(0),
            const S k4 = S(0))
        {
            p[0] = fx;
            p[1] = fy;
            p[2] = cx;
            p[3] = cy;
            p[4] = k1;
            p[5] = k2;
            p[6] = k3;
            p[7] = k4;
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
            
            const S r = sqrt(xp * xp + yp * yp);
            const S t = atan(r);
            const S t2 = t * t;
            const S t4 = t2 * t2;
            const S t6 = t4 * t2;
            const S t8 = t4 * t4;

            const S poly = (1 + p[K1] * t2 + p[K2] * t4 + p[K3] * t6 + p[K4] * t8);
            const S ts = t * poly;

            const S invr = (r < std::numeric_limits<S>::epsilon()) ? S(1) : S(1) / r;
            const S s = ts * invr;

            const S sxp = s * xp;
            const S syp = s * yp;

            u = sxp * p[FX] + p[CX];
            v = syp * p[FY] + p[CY];

            if (Jp)
            {
                S* Ju = Jp;
                S* Jv = Jp + NUM_PARAMS;

                Ju[FX] = sxp;
                Jv[FX] = S(0);

                Ju[FY] = S(0);
                Jv[FY] = syp;

                Ju[CX] = S(1);
                Jv[CX] = S(0);

                Ju[CY] = S(0);
                Jv[CY] = S(1);

                const S txpfx_r = t * invr * xp * p[FX];
                const S typfy_r = t * invr * yp * p[FY];

                Ju[K1] = t2 * txpfx_r;
                Jv[K1] = t2 * typfy_r;

                Ju[K2] = t4 * txpfx_r;
                Jv[K2] = t4 * typfy_r;

                Ju[K3] = t6 * txpfx_r;
                Jv[K3] = t6 * typfy_r;

                Ju[K4] = t8 * txpfx_r;
                Jv[K4] = t8 * typfy_r;
            }

            if (Jxyz)
            {
                S* Ju = Jxyz;
                S* Jv = Jxyz + 3;

                const S invr2 = invr * invr;
                const S r2p1 = r * r + 1;
                const S invr2p1 = (r2p1 < std::numeric_limits<S>::epsilon()) ? S(1) : S(1) / r2p1;

                const S ds_dr =
                    ((2 * p[K1] * t2) + (4 * p[K2] * t4) + (6 * p[K3] * t6) + (8 * p[K4] * t8)) * invr * invr2p1 -
                    ts * invr2 +
                    poly * invr * invr2p1;
                
                const S dr_dx = xp * invz * invr;
                const S dr_dy = yp * invz * invr;
                const S dr_dz = -r * invz;

                const S ds_dx = ds_dr * dr_dx;
                const S ds_dy = ds_dr * dr_dy;
                const S ds_dz = ds_dr * dr_dz;

                const S xpfx = xp * p[FX];
                const S ypfy = yp * p[FY];

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
            yz = (v - p[CY]) / p[FY];
            xz = (u - p[CX]) / p[FX];

            return IterativeUnproject(u, v, xz, yz, 10);
        }

        inline S& Param(int i) override { return p[i]; }
        inline const S& Param(int i) const override  { return p[i]; }
        inline int NumParams() const override { return NUM_PARAMS; }

        static std::string ModelName() { return "equidistant"; }
        inline std::string Name() const override { return ModelName(); }

        std::unique_ptr<CameraModel> Clone() const override
        {
            std::unique_ptr<CameraModel<S>> clone(new EquidistantModel<S>());
            clone->LoadModel(*this);
            return clone;
        }

        std::array<S, NUM_PARAMS> p;
    };
}
