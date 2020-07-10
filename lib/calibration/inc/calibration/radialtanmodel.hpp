#pragma once

#include <array>
#include <calibration/cameramodel.hpp>

namespace sight
{
    template <typename S>
    class RadialTanModel : public CameraModel<S>
    {
    public:
        enum
        {
            FX = 0,
            FY,
            CX,
            CY,
            P1,
            P2,
            K1,
            K2,
            K3,
            K4,
            K5,
            K6,
            NUM_PARAMS
        };

        RadialTanModel(
            const S fx = S(1),
            const S fy = S(1),
            const S cx = S(0),
            const S cy = S(0),
            const S p1 = S(0),
            const S p2 = S(0),
            const S k1 = S(0),
            const S k2 = S(0),
            const S k3 = S(0),
            const S k4 = S(0),
            const S k5 = S(0),
            const S k6 = S(0))
        {
            p[FX] = fx;
            p[FY] = fy;
            p[CX] = cx;
            p[CY] = cy;
            p[P1] = p1;
            p[P2] = p2;
            p[K1] = k1;
            p[K2] = k2;
            p[K3] = k3;
            p[K4] = k4;
            p[K5] = k5;
            p[K6] = k6;
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

            const S xp2 = xp * xp;
            const S yp2 = yp * yp;
            const S xpyp = xp * yp;

            const S r2 = xp2 + yp2;
            const S r4 = r2 * r2;
            const S r6 = r4 * r2;

            const S num = S(1) + (p[K1] * r2) + (p[K2] * r4) + (p[K3] * r6);
            const S den = S(1) + (p[K4] * r2) + (p[K5] * r4) + (p[K6] * r6);

            const S deninv = (abs(den) > std::numeric_limits<S>::epsilon()) ? S(1) / den : S(1);
            const S s = num * deninv;
            
            const S xpp = xp * s + (2 * p[P1] * xpyp) + (p[P2] * (r2 + 2 * xp2));
            const S ypp = yp * s + (p[P1] * (r2 + 2 * yp2)) + (2 * p[P2] * xpyp);

            u = p[FX] * xpp + p[CX];
            v = p[FY] * ypp + p[CY];

            if (Jp)
            {
                S* Ju = Jp;
                S* Jv = Jp + NUM_PARAMS;

                Ju[FX] = xpp;
                Jv[FX] = S(0);

                Ju[FY] = S(0);
                Jv[FY] = ypp;

                Ju[CX] = S(1);
                Jv[CX] = S(0);

                Ju[CY] = S(0);
                Jv[CY] = S(1);

                Ju[P1] = p[FX] * (2 * xpyp);
                Jv[P1] = p[FY] * (r2 + 2 * yp2);

                Ju[P2] = p[FX] * (r2 + 2 * xp2);
                Jv[P2] = p[FY] * (2 * xpyp);

                const S fxxp = p[FX] * xp;
                const S fyyp = p[FY] * yp;

                const S r2_den = r2 * deninv;
                Ju[K1] = fxxp * r2_den;
                Jv[K1] = fyyp * r2_den;

                const S r4_den = r4 * deninv;
                Ju[K2] = fxxp * r4_den;
                Jv[K2] = fyyp * r4_den;

                const S r6_den = r6 * deninv;
                Ju[K3] = fxxp * r6_den;
                Jv[K3] = fyyp * r6_den;

                const S s_den = s * deninv;

                const S sr2_den = s_den * r2;
                Ju[K4] = fxxp * -sr2_den;
                Jv[K4] = fyyp * -sr2_den;

                const S sr4_den = s_den * r4;
                Ju[K5] = fxxp * -sr4_den;
                Jv[K5] = fyyp * -sr4_den;

                const S sr6_den = s_den * r6;
                Ju[K6] = fxxp * -sr6_den;
                Jv[K6] = fyyp * -sr6_den;
            }

            if (Jxyz)
            {
                S* Ju = Jxyz;
                S* Jv = Jxyz + 3;

                const S r = sqrt(r2);
                const S invr = (r > std::numeric_limits<S>::epsilon()) ? S(1) / r : S(1);

                const S dnum_dr = r * ((S(2) * p[K1]) + (S(4) * p[K2] * r2) + (S(6) * p[K3] * r4));
                const S dden_dr = r * ((S(2) * p[K4]) + (S(4) * p[K5] * r2) + (S(6) * p[K6] * r4));
                const S ds_dr = deninv * (dnum_dr - dden_dr * s);

                const S invrz = invz * invr;
                const S dr_dx = xp * invrz;
                const S dr_dy = yp * invrz;
                const S dr_dz = -r * invz;

                const S ds_dx = ds_dr * dr_dx;
                const S ds_dy = ds_dr * dr_dy;
                const S ds_dz = ds_dr * dr_dz;

                const S xp_z = xp * invz;
                const S yp_z = yp * invz;

                Ju[0] = p[FX] * ((invz * s + xp * ds_dx) + (S(2) * p[P1] * yp_z) + (p[P2] * S(6) * xp_z));
                Jv[0] = p[FY] * ((yp * ds_dx) + (S(2) * p[P1] * xp_z) + (S(2) * p[P2] * yp_z));

                Ju[1] = p[FX] * ((xp * ds_dy) + (S(2) * p[P1] * xp_z) + (S(2) * p[P2] * yp_z));
                Jv[1] = p[FY] * ((invz * s + yp * ds_dy) + (S(6) * p[P1] * yp_z) + (p[P2] * S(2) * xp_z));

                Ju[2] = p[FX] * ((-xp * invz * s + xp * ds_dz) + (S(-4) * p[P1] * xpyp * invz) + (S(-2) * p[P2] * invz * (S(3) * xp2 + yp2)));
                Jv[2] = p[FY] * ((-yp * invz * s + yp * ds_dz) + (S(-2) * p[P1] * invz * (xp2 + S(3) * yp2)) + (S(-4) * p[P2] * xpyp * invz));
            }

            return true;
        }

        bool Unproject(S u, S v, S& xz, S& yz) const override
        {
            yz = (v - p[CY]) / p[FY];
            xz = (u - p[CX]) / p[FX];

            return IterativeUnproject(u, v, xz, yz, 100);
        }

        inline S& Param(int i) override { return p[i]; }
        inline const S& Param(int i) const override  { return p[i]; }
        inline int NumParams() const override { return NUM_PARAMS; }

        inline std::string Name() const override { return "radialtan"; }

        CameraModel* Clone() const override
        {
            RadialTanModel* clone = new RadialTanModel<S>(p[FX]);
            std::copy(p.begin(), p.end(), clone->p.begin());
            return clone;
        }

        std::array<S, NUM_PARAMS> p;
    };
}
