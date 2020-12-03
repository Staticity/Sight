#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

#include <linear/vec.hpp>
#include <geometric/homography_ops.hpp>
#include <geometric/rectify.hpp>
#include <image/image_ops.hpp>

// An implementation of Marc Pollefeys' polar rectification
// https://people.inf.ethz.ch/pomarc/pubs/PollefeysICCV99.pdf
//

namespace sight
{

    namespace
    {
        // Helper functions placed in an anonymous namespace
        enum Border
        {
            TOP = 0,
            BOTTOM = 1,
            LEFT = 2,
            RIGHT = 3
        };

        // Returns an integer in the range [1, 9] representing
        // where the epipole lies around [0, 1] x [0, 1].
        //
        //    1 | 2 | 3
        //  ----a---b----
        //    4 | 5 | 6
        //  ----c---d----
        //    7 | 8 | 9
        template <typename S>
        int ComputeEpipolarRegion(const Vec2<S>& e, int w, int h)
        {
            const int xi = (e(0) < S(0)) ? 0 : (e(0) < S(w)) ? 1 : 2;
            const int yi = (e(1) < S(0)) ? 0 : (e(1) < S(h)) ? 1 : 2;
            return yi * 3 + xi + 1;
        }

        template <typename S>
        bool EpipoleInsideImage(const Vec2<S>& e, int w, int h)
        {
            return ComputeEpipolarRegion(e, w, h) == 5;
        }

        // Stores the extremal corners based on the region index provided.
        // Returns false if the epipole is inside the image.
        template <typename S>
        bool GetExtremalCorners(int region, int w, int h, Vec2<S>& c0, Vec2<S>& c1)
        {
            const Vec2<S> a = {S(0), S(0)};
            const Vec2<S> b = {S(w), S(0)};
            const Vec2<S> c = {S(0), S(h)};
            const Vec2<S> d = {S(w), S(h)};

            switch (region)
            {
            case 1:
                c0 = b;
                c1 = c;
                return true;
            case 2:
                c0 = b;
                c1 = a;
                return true;
            case 3:
                c0 = d;
                c1 = a;
                return true;
            case 4:
                c0 = a;
                c1 = c;
                return true;
            case 6:
                c0 = d;
                c1 = b;
                return true;
            case 7:
                c0 = a;
                c1 = d;
                return true;
            case 8:
                c0 = c;
                c1 = d;
                return true;
            case 9:
                c0 = c;
                c1 = b;
                return true;
            default:
                return false;
            }
        }

        template <typename S>
        void GetPolarMinMaxAngles(
            const Vec2<S>& e,
            int w,
            int h,
            S& theta0,
            S& theta1)
        {
            // TODO: Something is up with the direction of the angles.. Maybe becaues image space has -y???
            Vec2<S> c0, c1;
            const int region = ComputeEpipolarRegion(e, w, h);
            const bool outsideImage = GetExtremalCorners(region, w, h, c0, c1);
            if (outsideImage)
            {
                const Vec2<S> d0 = c0 - e;
                const Vec2<S> d1 = c1 - e;

                theta0 = atan2(d0(1), d0(0));
                theta1 = atan2(d1(1), d1(0));

                // Bring angles into the range [0, 2pi]
                if (theta0 < S(0))
                {
                    theta0 += 2 * M_PI;
                }

                if (theta1 < S(0))
                {
                    theta1 += 2 * M_PI;
                }
            }
            else
            {
                theta0 = S(0);
                theta1 = S(2 * M_PI);
            }
        }

        template <typename S>
        void GetImageBorderLines(Vec3<S> lines[4], int w, int h)
        {
            // Should these be in image space...?
            lines[Border::TOP] = Vec3<S>({S(0), S(0), S(1)}).Cross({S(w), S(0), S(1)});
            lines[Border::BOTTOM] = Vec3<S>({S(0), S(h), S(1)}).Cross({S(w), S(h), S(1)});

            lines[Border::LEFT] = Vec3<S>({S(0), S(0), S(1)}).Cross({S(0), S(h), S(1)});
            lines[Border::RIGHT] = Vec3<S>({S(w), S(0), S(1)}).Cross({S(w), S(h), S(1)});
        }

        template <typename S>
        void GetPolarMinMaxRho(
            const Vec2<S>& e,
            const Vec2<S>& v,
            const Vec3<S> borders[4],
            int w,
            int h,
            S& minR,
            S& maxR)
        {
            minR = std::numeric_limits<S>::max();
            maxR = std::numeric_limits<S>::min();

            // If we're in the image, we need to include the epipole.
            if (EpipoleInsideImage(e, w, h))
            {
                minR = S(0);
                maxR = S(0);
            }

            // Check for intersection among the borders and only accept
            // those with positive values of theta.
            const Vec3<S> eh = ToHomogeneous(e);
            const Vec3<S> e_off = { eh(0) + v(0), eh(1) + v(1), eh(2) };

            // We normalize this to keep the scale uniform across
            // all calls to this function.
            const Vec2<S> vn = v.Normalized();
            
            const Vec3<S> line = eh.Cross(e_off);
            for (int i = 0; i < 4; ++i)
            {
                // Intersection of two lines is the cross product
                Vec3<S> pth = borders[i].Cross(line);

                // Point at infinity check -- aka, lines are parallel
                if (abs(pth(2)) < std::numeric_limits<S>::epsilon())
                {
                    continue;
                }

                const Vec2<S> pt = FromHomogeneous(pth);

                // Exclude points outside the rectangle
                if (pt(0) < S(0) || pt(0) > S(w) ||
                    pt(1) < S(0) || pt(1) > S(h))
                {
                    continue;
                }

                // We have a ray of the form:
                //
                //     pt = e + t * v
                //
                // where v = (cos(theta), sin(theta))
                //
                // We can solve for t by finding the t which
                // minimizes the squared distance from the point.
                //
                // So minimize (e + t*v - pt)^T(e + t*v - pt)
                //
                // Taking the derivative w.r.t to t and setting to 0
                // results in:
                //
                //     v^T * (e + t*v - pt) = 0
                //
                // which implies:
                //
                //    t = v^T * (pt - e) / (v^T * v)
                
                const S t = vn.Dot(pt - e) / vn.SquaredNorm();

                // Ignore points going in the opposite direction
                if (t < S(0))
                {
                    continue;
                }

                minR = std::min(t, minR);
                maxR = std::max(t, maxR);
            }
        }

        template <typename S>
        Eigen::Matrix3<S> Skew(const Vec3<S>& v)
        {
            Eigen::Matrix3<S> skew;
            skew <<
                S(0), -v(2), v(1),
                v(2), S(0), -v(0),
                -v(1), v(0), S(0);
            return skew;
        }

        template <typename S>
        bool InBetweenAngles(const S t0, const S t1, S theta)
        {
            while (theta < t0)
            {
                theta += 2 * M_PI;
            }

            while (theta > t1)
            {
                theta -= 2 * M_PI;
            }

            return theta >= t0 && theta <= t1;
        }
    } // namespace

    template <typename S>
    class InversePolarWarp : public IInverseWarp<S>
    {
    public:

        InversePolarWarp()
            : uMax(std::numeric_limits<int>::min())
            , vMax(std::numeric_limits<int>::min())
        {}

        bool InterpolateUV(S x, S y, S& u, S& v) const
        {
            if (x < S(0) || x + 1 >= uvFromRhoTheta.w ||
                y < S(0) || y + 1 >= uvFromRhoTheta.h)
            {
                return false;
            }

            u = BilinearInterpolate<S, S, S>(uvFromRhoTheta, x, y, 0);
            v = BilinearInterpolate<S, S, S>(uvFromRhoTheta, x, y, 1);

            return true;
        }
        
        void GetUV(int x, int y, S& u, S& v) const
        {
            if (x < S(0) || x + 1 > uvFromRhoTheta.w ||
                y < S(0) || y + 1 > uvFromRhoTheta.h)
            {
                u = S(0);
                v = S(0);
                return;
            }
            
            const S* ptr = uvFromRhoTheta.at(y, x, 0);
            u = ptr[0];
            v = ptr[1];
        }

        virtual void GetUVBounds(int& w, int& h) const
        {
            w = uMax;
            h = vMax;
        }

        int uMax;
        int vMax;
        Image<S> uvFromRhoTheta;
    };

    template <typename T, typename S>
    void PolarRectifyImages(
        const Image<T>& im0,
        const CameraModel<S>& pinhole0,
        const Image<T>& im1,
        const CameraModel<S>& pinhole1,
        const SE3<S>& cam1FromCam0,
        Image<T>& out0,
        InversePolarWarp<S>& invWarp0,
        Image<T>& out1,
        InversePolarWarp<S>& invWarp1)
    {
        // assert(pinhole0.Name() == "pinhole")
        // assert(pinhole1.Name() == "pinhole")
        Vec3<S> e0h, e1h;
        Eigen::Matrix3<S> F;

        // All of our camera models use normalized pixel space. So we'll undo
        // that to keep everything simpler in pixel space in this function.
        auto K0 = pinhole0.Clone();
        auto K1 = pinhole0.Clone();
        {
            K0->Param(PinholeModel<S>::FX) *= im0.w;
            K0->Param(PinholeModel<S>::FY) *= im0.h;
            K0->Param(PinholeModel<S>::CX) *= im0.w;
            K0->Param(PinholeModel<S>::CY) *= im0.h;

            K1->Param(PinholeModel<S>::FX) *= im1.w;
            K1->Param(PinholeModel<S>::FY) *= im1.h;
            K1->Param(PinholeModel<S>::CX) *= im1.w;
            K1->Param(PinholeModel<S>::CY) *= im1.h;
        }
        ComputeEpipolesFromLinearCalibration(*K0, *K1, cam1FromCam0, F, e0h, e1h);

        const Vec2<S> e0 = FromHomogeneous(e0h);
        const Vec2<S> e1 = FromHomogeneous(e1h);

        S tMin0, tMax0;
        GetPolarMinMaxAngles<S>(e0, im0.w, im0.h, tMin0, tMax0);

        S tMin1, tMax1;
        GetPolarMinMaxAngles<S>(e1, im1.w, im1.h, tMin1, tMax1);

        S thetaMin = std::max(tMin0, tMin1);
        S thetaMax = std::min(tMax0, tMax1);

        // Ensure thetaMax > thetaMin (should only run once)
        while (thetaMax < thetaMin)
        {
            thetaMax += 2 * M_PI;
        }

        S theta_left = thetaMin;

        Vec3<S> border0[4], border1[4];
        GetImageBorderLines(border0, im0.w, im0.h);
        GetImageBorderLines(border1, im1.w, im1.h);

        // This is not the correct dTheta. Technically the
        // maximum nThetas should be 2 * (w + h), but we leave
        // this as is for performance reasons...
        const int nThetas = (im0.w + im0.h);
        const S dTheta = (thetaMax - thetaMin) / (nThetas - 1);

        std::vector<std::pair<S, S>> rho_starts;
        std::vector<std::pair<Vec2<S>, Vec2<S>>> v_pairs;
        rho_starts.reserve(nThetas);
        v_pairs.reserve(nThetas);

        S rhoMaxRange = std::numeric_limits<S>::min();

        for (int i = 0; i < nThetas; ++i)
        {
            S rMin0, rMax0;
            const Vec2<S> v0 = { cos(theta_left), sin(theta_left) };
            GetPolarMinMaxRho<S>(e0, v0, border0, im0.w, im0.h, rMin0, rMax0);

            // Compute theta in the right-image.
            //
            // Find the epipolar line corresponding to test_theta in the left & right image
            const S r = rMax0;
            const Vec3<S> off = { r * v0(0), r * v0(1), S(0) };
            const Eigen::Vector3d eLine0 = ToEigen(e0h.Cross(e0h + off));
            const Eigen::Vector3d eLine1 = F * Skew(e0h) * eLine0;

            // Ensure this is a right-epipolar line (passes through right-epipole)
            assert(abs(eLine1.transpose() * ToEigen(e1h)) < std::numeric_limits<S>::epsilon());

            const S c = eLine1(2);
            const S a = eLine1(0) / c;
            const S b = eLine1(1) / c;

            // TODO: Not sure if this calculation for `theta_right` is always accurate..
            S rMin1, rMax1;
            const S theta_right = atan2(-a, b);
            const Vec2<S> v1 = { cos(theta_right), sin(theta_right) };
            GetPolarMinMaxRho<S>(e1, v1, border1, im1.w, im1.h, rMin1, rMax1);

            // Store the bounds of rho
            rho_starts.push_back({rMin0, rMin1});
            rhoMaxRange = std::max({rMax0 - rMin0, rMax1 - rMin1, rhoMaxRange});

            // Store the direction of each epipolar line
            v_pairs.push_back({v0, v1});
            
            // TODO: Compute proper dTheta
            theta_left += dTheta;
        }

        const int nRhos = int(ceil(rhoMaxRange));

        out0 = Image<T>(nRhos, nThetas, im0.c);
        out1 = Image<T>(nRhos, nThetas, im1.c);
        invWarp0.uvFromRhoTheta = Image<S>(nRhos, nThetas, 2);
        invWarp1.uvFromRhoTheta = Image<S>(nRhos, nThetas, 2);

        // Set inverse warp's valid bounds
        invWarp0.uMax = im0.w;
        invWarp0.vMax = im0.h;

        invWarp1.uMax = im1.w;
        invWarp1.vMax = im1.h;

        const S thetaRange = (thetaMax - thetaMin);
        for (int i = 0; i < nThetas; ++i)
        {
            const auto& v0 = v_pairs[i].first;
            const auto& v1 = v_pairs[i].second;

            const auto& rStart = rho_starts[i];
            for (int j = 0; j < nRhos; ++j)
            {
                const S rho0 = rStart.first + S(j);
                const S rho1 = rStart.second + S(j);
                auto uv0 = e0 + v0 * rho0;
                auto uv1 = e1 + v1 * rho1;

                invWarp0.uvFromRhoTheta(i, j, 0) = uv0(0);
                invWarp0.uvFromRhoTheta(i, j, 1) = uv0(1);

                invWarp1.uvFromRhoTheta(i, j, 0) = uv1(0);
                invWarp1.uvFromRhoTheta(i, j, 1) = uv1(1);

                // Can we bilinear interpolate this pixel?
                if (uv0[0] > S(0) && uv0[0] + S(1) < im0.w &&
                    uv0[1] > S(0) && uv0[1] + S(1) < im0.h)
                {
                    for (int ch = 0; ch < im0.c; ++ch)
                    {
                        out0(i, j, ch) = BilinearInterpolate<T, T, S>(im0, uv0[0], uv0[1], ch);
                    }
                }

                if (uv1[0] > S(0) && uv1[0] + S(1) < im1.w &&
                    uv1[1] > S(0) && uv1[1] + S(1) < im1.h)
                {
                    for (int ch = 0; ch < im1.c; ++ch)
                    {
                        out1(i, j, ch) = BilinearInterpolate<T, T, S>(im1, uv1[0], uv1[1], ch);
                    }
                }
            }
        }
    }
}
