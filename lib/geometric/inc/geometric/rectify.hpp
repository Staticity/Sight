#pragma once

#include <random>
#include <time.h>

#include <Eigen/Eigen>

#include <calibration/pinholemodel.hpp>
#include <calibration/device.hpp>
#include <estimation/quadratic_spline.hpp>
#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <image/image_gen.hpp>
#include <geometric/homography_ops.hpp>
#include <geometric/triangulation.hpp>

namespace sight
{
    template <typename S>
    class IInverseWarp
    {
    public:

        virtual bool InterpolateUV(S x, S y, S& u, S& v) const = 0;
        
        virtual void GetUV(int x, int y, S& u, S& v) const = 0;

        virtual void GetUVBounds(int& w, int& h) const = 0;
    };

    template <typename S>
    void DecomposeFundamentalIntoSM(
        const Eigen::Matrix3<S>& F,
        // Eigen::Matrix3<S>& S,
        // This is the skew symmetric matrix of the left epipole.
        // We ignore it for now.
        Eigen::Matrix3<S>& M)
    {
        const Eigen::JacobiSVD svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const auto sv = svd.singularValues();

        Eigen::Matrix3<S> Lp;
        Lp <<
            sv(0), S(0), S(0),
            S(0), sv(1), S(0),
            S(0), S(0), (sv(0) + sv(1)) / S(2);
        
        Eigen::Matrix3<S> W;
        W <<
            S(0), S(-1), S(0),
            S(1), S(0), S(0),
            S(0), S(0), S(1);
        
        M = svd.matrixU() * W * Lp * svd.matrixV().transpose();
    }

    template <typename S>
    void ComputeEpipolesFromEssentialOrFundamental(
        const Eigen::Matrix3<S>& F,
        S e0[3],
        S e1[3])
    {
        const Eigen::JacobiSVD svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::Vector3<S> epi0 = svd.matrixV().col(2);
        const Eigen::Vector3<S> epi1 = svd.matrixU().col(2);

        e0[0] = epi0[0];
        e0[1] = epi0[1];
        e0[2] = epi0[2];

        e1[0] = epi1[0];
        e1[1] = epi1[1];
        e1[2] = epi1[2];
    }

    template <typename S>
    Eigen::Matrix3<S> ComputeEssentialFromRelativePose(
        const SE3<S>& cam1FromCam0)
    {
        const auto& Rot = cam1FromCam0.R;
        const auto& t = cam1FromCam0.t;
        Eigen::Matrix3<S> R;
        R <<
            Rot[0], Rot[1], Rot[2],
            Rot[3], Rot[4], Rot[5],
            Rot[6], Rot[7], Rot[8];

        Eigen::Matrix3<S> Tx;
        Tx <<
             S(0), -t[2],  t[1],
             t[2],  S(0), -t[0],
            -t[1],  t[0],  S(0);

        return Tx * R;
    }

    template <typename S>
    Eigen::Matrix3<S> ComputeFundamentalFromCalibration(
        const CameraModel<S>& cam0,
        const CameraModel<S>& cam1,
        const SE3<S>& cam1FromCam0)
    {
        assert(cam0.Name() == "pinhole");
        assert(cam1.Name() == "pinhole");

        // Linear calibrations are 3x3 matrices of the form:
        //
        // [fx 0 cx]
        // [0 fy cy]
        // [0  0  1]
        //
        // With inverses:
        //
        // [1/fx     0  -cx/fx]
        // [   0  1/fy  -cy/fy]
        // [   0     0       1]

        Eigen::Matrix3<S> K0_inv;

        // Compute K0^-1
        {
            const S fx = cam0.Param(PinholeModel<S>::FX);
            const S fy = cam0.Param(PinholeModel<S>::FY);
            const S cx = cam0.Param(PinholeModel<S>::CX);
            const S cy = cam0.Param(PinholeModel<S>::CY);

            K0_inv <<
                S(1) / fx,      S(0), -cx / fx,
                S(0)     , S(1) / fy, -cy / fy,
                S(0)     ,      S(0),     S(1);
        }

        Eigen::Matrix3<S> K1_inv;
        // Compute K1^-1
        {
            const S fx = cam1.Param(PinholeModel<S>::FX);
            const S fy = cam1.Param(PinholeModel<S>::FY);
            const S cx = cam1.Param(PinholeModel<S>::CX);
            const S cy = cam1.Param(PinholeModel<S>::CY);

            K1_inv <<
                S(1) / fx,      S(0), -cx / fx,
                S(0)     , S(1) / fy, -cy / fy,
                S(0)     ,      S(0),     S(1);
        }

        const Eigen::Matrix3<S> E = ComputeEssentialFromRelativePose(cam1FromCam0);
        return K1_inv.transpose() * E * K0_inv;
    }

    template <typename S>
    void ComputeEpipolesFromLinearCalibration(
        const CameraModel<S>& cam0,
        const CameraModel<S>& cam1,
        const SE3<S>& cam1FromCam0,
        Eigen::Matrix3<S>& F,
        Vec3<S>& e0,
        Vec3<S>& e1)
    {
        assert(cam0.Name() == "pinhole");
        assert(cam1.Name() == "pinhole");

        F = ComputeFundamentalFromCalibration(cam0, cam1, cam1FromCam0);
        ComputeEpipolesFromEssentialOrFundamental(F, &e0[0], &e1[0]);
    }

    template <typename S>
    void ComputeEpipolesFromCalibration(
        const CameraModel<S>& cam0,
        const CameraModel<S>& cam1,
        const SE3<S>& cam1FromCam0,
        Vec3<S>& e0,
        Vec3<S>& e1)
    {
        const Eigen::Matrix3<S> E = ComputeEssentialFromRelativePose(cam1FromCam0);
        
        Vec3<S> me0, me1;
        ComputeEpipolesFromEssentialOrFundamental(E, &me0[0], &me1[0]);

        // Project the epipolar metric rays back into the images
        cam0.Project(me0[0], me0[1], me0[2], e0[0], e0[1]);
        e0[2] = S(1);

        cam1.Project(me1[0], me1[1], me1[2], e1[0], e1[1]);
        e1[2] = S(1);
    }

    template <typename S>
    void ComputeEpipolesFromCalibration(
        const Camera<S>& cam0,
        const Camera<S>& cam1,
        Vec3<S>& e0,
        Vec3<S>& e1)
    {
        ComputeEpipolesFromCalibration(*cam0.model, *cam1.model, cam1.Rt * cam0.Rt.Inverse());
    }

    template <typename S>
    PinholeModel<S> FindOptimalLinearCalibration(
        const CameraModel<S>& cam,
        const int currentWidth,
        const int currentHeight,
        const int desiredWidth,
        const int desiredHeight,
        const S xRadius = std::numeric_limits<S>::max(),
        const S yRadius = std::numeric_limits<S>::max())
    {
        // Let's unproject all of the boundary pixels to
        // find the bounding rectangle of all the coordinates.
        const S invW = S(1) / currentWidth;
        const S invH = S(1) / currentHeight;
        S minX = std::numeric_limits<S>::max();
        S minY = std::numeric_limits<S>::max();
        S maxX = std::numeric_limits<S>::min();
        S maxY = std::numeric_limits<S>::min();

        // We choose to unproject the entirety of the image
        // insetad of just the boundary, since some camera
        // models may not be able to unproject the boundaries
        // due to extreme distortion.
        //
        // We will slightly optimize this process by scanning
        // each scanline left-to-right until success, and vice-versa.
        for (int i = 0; i < currentHeight; ++i)
        {
            S xz, yz;

            int l;
            for (l = 0; l < currentWidth; ++l)
            {
                if (cam.Unproject((S(l) + S(0.5)) * invW, (S(i) + S(0.5)) * invH, xz, yz))
                {
                    minX = std::min(xz, minX);
                    minY = std::min(yz, minY);
                    maxX = std::max(xz, maxX);
                    maxY = std::max(yz, maxY);
                    break;
                }
            }

            for (int r = currentWidth - 1; r > l; --r)
            {
                if (cam.Unproject((S(r) + S(0.5)) * invW, (S(i) + S(0.5)) * invH, xz, yz))
                {
                    minX = std::min(xz, minX);
                    minY = std::min(yz, minY);
                    maxX = std::max(xz, maxX);
                    maxY = std::max(yz, maxY);
                    break;
                }
            }
        }

        // Clamp the bounding box to the maxRadius defined
        minX = std::max(minX, -xRadius);
        minY = std::max(minY, -xRadius);
        maxX = std::min(maxX, xRadius);
        maxY = std::min(maxY, xRadius);

        // Now that we have the bounds of the undistortion, we can
        // solve for a camera calibration which projects the bounding
        // rectangle into the bounds of our image.
        //
        // So find K such that K * {minX, minY, 1} = (0.5 / w, 0.5 / h)
        // and so forth.
        //
        // K has 5 parameters:
        //
        //     fx, fy, cx, cy, and skew
        //
        // fx and fy are simply the ratio between the bounding rectangle
        // and the new image's bounds.
        //
        // Then we have:
        //
        //     fx * minX + cx = 0.5 / w
        //     fy * minY + cy = 0.5 / h
        //
        // So, 

        const S dx = maxX - minX;
        const S dy = maxY - minY;

        const S fx = desiredWidth / dx;
        const S fy = desiredHeight / dy;
        const S cx = S(0.5) * invW - fx * minX;
        const S cy = S(0.5) * invH - fy * minY;

        const S destInvW = S(1) / desiredWidth;
        const S destInvH = S(1) / desiredHeight;

        PinholeModel<S> outCam;
        outCam.Param(PinholeModel<S>::FX) = fx * destInvW;
        outCam.Param(PinholeModel<S>::FY) = fy * destInvH;
        outCam.Param(PinholeModel<S>::CX) = cx * destInvW;
        outCam.Param(PinholeModel<S>::CY) = cy * destInvH;
        outCam.Param(PinholeModel<S>::SKEW) = S(0);
        return outCam;
    }

    template <typename T, typename S>
    void UndistortImage(
        const Image<T>& im,
        const CameraModel<S>& cam,
        const PinholeModel<S>& newCam,
        Image<T>& outIm)
    {
        if (outIm.w == 0 || outIm.h == 0)
        {
            outIm = Image<T>(im.w, im.h, im.c);
        }

        // Now, we must inverse sample the pixels in the original
        // image into the new image.
        const S invW = S(1) / outIm.w;
        const S invH = S(1) / outIm.h;
        auto view = PadView(im, 1);
        for (int i = 0; i < outIm.h; ++i)
        {
            T* row = outIm.row(i);
            for (int j = 0, k = 0; j < outIm.w; ++j)
            {
                bool ok = true;
                // Map back into metric space
                S xz, yz;
                ok &= newCam.Unproject(
                    (S(j) + S(.5)) * invW,
                    (S(i) + S(.5)) * invH,
                    xz,
                    yz);

                // Project into the original image
                S u, v;
                ok &= cam.Project(xz, yz, S(1), u, v, nullptr, nullptr);
                u = (u * im.w) - S(.5);
                v = (v * im.h) - S(.5);

                ok &= u >= 0 && u <= view.w;
                ok &= v >= 0 && v <= view.h;

                // Bilinear sample this location for each location
                if (ok)
                {
                    for (int ch = 0; ch < outIm.c; ++ch, ++k)
                    {
                        row[k] = BilinearInterpolate<T, T, S>(view, u, v, ch);
                    }
                }
                else
                {
                    for (int ch = 0; ch < outIm.c; ++ch, ++k)
                    {
                        row[k] = T(0);
                    }
                }
            }
        }
    }

    template <typename T, typename S>
    void UndistortImage(
        const Image<T>& im,
        const CameraModel<S>& cam,
        Image<T>& outIm,
        PinholeModel<S>& outCam,
        const S xRadius = std::numeric_limits<S>::max(),
        const S yRadius = std::numeric_limits<S>::max())
    {
        if (outIm.w == 0 || outIm.h == 0)
        {
            outIm = Image<T>(im.w, im.h, im.c);
        }

        outCam = FindOptimalLinearCalibration(cam, im.w, im.h, outIm.w, outIm.h, xRadius, yRadius);
        UndistortImage(im, cam, outCam, outIm);
    }

    template <typename T, typename S>
    void RemapPerspectiveImage(
        const Image<T>& im,
        const Eigen::Matrix3<S>& H,
        Image<T>& out)
    {
        // Determine the bounding box by warping the extremities
        S minX = std::numeric_limits<S>::max();
        S minY = std::numeric_limits<S>::max();
        S maxX = std::numeric_limits<S>::min();
        S maxY = std::numeric_limits<S>::min();

        // Iterate over the corners to find the extrema
        for (int i = 0; i < 4; ++i)
        {
            const int x = (i % 2) ? 0 : im.w;
            const int y = (i / 2) ? 0 : im.h;

            const Eigen::Vector2<S> uv = ApplyHomography(H, Eigen::Vector2<S>(x, y));
            minX = std::min(uv[0], minX);
            minY = std::min(uv[1], minY);
            maxX = std::max(uv[0], maxX);
            maxY = std::max(uv[1], maxY);
        }

        // Inverse warp into the original image to sample
        const Eigen::Matrix3<S> Hinv = H.inverse();

        const int w = int(ceil(maxX - minX));
        const int h = int(ceil(maxY - minY));
        out = Image<T>(w, h, im.c);

        for (int i = 0; i < out.h; ++i)
        {
            T* dstRow = out.row(i);
            for (int j = 0, k = 0; j < out.w; ++j)
            {
                const Eigen::Vector2<S> uvDst = Eigen::Vector2d(S(j), S(i));
                const Eigen::Vector2<S> uvSrc = ApplyHomography(Hinv, uvDst);

                // Skip pixels which are out of the image's bounds
                if (uvSrc[0] < S(0) || uvSrc[0] >= out.w - 1 ||
                    uvSrc[1] < S(0) || uvSrc[1] >= out.h - 1)
                {
                    continue;
                }

                for (int ch = 0; ch < out.c; ++ch, ++k)
                {
                    dstRow[k] = BilinearInterpolate<T, T, S>(im, uvSrc[0], uvSrc[1], ch);
                }
            }
        }
    }

    template <typename T, typename S>
    Image<S> DisparityFromStereoRectified(
        const Image<T>& rectified0,
        const Image<T>& rectified1,
        const int windowRadius = 2,
        const float maxDisparityPerc = .01f,
        const bool refineSubpixel = true)
    {
        assert(rectified0.c == 1);
        assert(rectified1.c == 1);

        const auto& rec0 = rectified0;
        const auto& rec1 = rectified1;

        // Since the images are rectified, we only need to perform
        // searches along their horizontal epipolar lines.
        const int maxDisparity = int(round(rec0.w * maxDisparityPerc));
        const int disparityRadius = (maxDisparity + 1) / 2;

        const int wr = windowRadius;
        const int dr = disparityRadius;

        // For each pixel, form a patch around it and search along the
        // epipolar line to find a patch which matches most.
        Image<S> disparity(rec0.w, rec0.h);

        const int size = wr * 2 + 1;
        std::vector<S> ssds;
        ssds.reserve(size * size);

        for (int y = wr; y < rec0.h - wr; ++y)
        {
            const T* row0 = rec0.row(y);
            const T* row1 = rec1.row(y);
            for (int x = wr; x < rec0.w - wr; ++x)
            {
                // Compare this pixel's (x, y) patch against a set of
                // patches in the other image via SSD.
                const int start = std::max(0, x - disparityRadius);
                const int end = std::min(rec0.w - wr - 1, x + disparityRadius);

                S d = std::numeric_limits<S>::max();
                S minSsd = std::numeric_limits<S>::max();
                
                int bestIdx = -1;
                ssds.clear();

                const T* window0 = row0 + -wr * rec0.row_step + x;

                // For every patch center in image 1...
                for (int xp = start; xp <= end; ++xp)
                {
                    S ssd = S(0);

                    const T* win0 = row1 + -wr * rec1.row_step + xp;
                    const T* win1 = row1 + -wr * rec1.row_step + xp;

                    // For every pixel in that patch...
                    for (int i = -wr; i <= wr; ++i)
                    {
                        for (int j = -wr; j <= wr; ++j)
                        {
                            // Assuming there's more than 1 channel...
                            // const S err =
                            //     row0[i * rec0.row_step + (x + j) * rec0.col_step] -
                            //     row1[i * rec1.row_step + (xp + j) * rec1.col_step];

                            // const S err =
                            //     row0[i * rec0.row_step + (x + j)] -
                            //     row1[i * rec1.row_step + (xp + j)];

                            const S err = win0[j] - win1[j];
                            ssd += err * err;
                        }

                        win0 += rec0.row_step;
                        win1 += rec1.row_step;
                    }

                    if (refineSubpixel)
                    {
                        // spline.AddPoint(xp, ssd);
                        ssds.push_back(ssd);
                    }

                    if (ssd < minSsd)
                    {
                        minSsd = ssd;
                        d = xp;
                        bestIdx = xp - start;
                    }
                }

                if (refineSubpixel && bestIdx > 0 && bestIdx + 1 < ssds.size())
                {
                    const int i = bestIdx;
                    Quadratic<S> quad(d - 1, ssds[i - 1], d, ssds[i], d + 1, ssds[i + 1]);
                    quad.GetExtremum(d);
                }

                disparity(y, x) = S(d - x);
            }
        }

        return disparity;
    }

    template <typename S>
    Image<S> DepthFromDisparity(
        const Image<S>& disparity,
        const CameraModel<S>& pinhole0,
        const IInverseWarp<S>& invWarp0,
        const CameraModel<S>& pinhole1,
        const IInverseWarp<S>& invWarp1,
        const SE3<S>& cam1FromCam0)
    {
        assert(pinhole0.Name() == "pinhole");
        assert(pinhole1.Name() == "pinhole");

        const SE3<S> cam0FromCam1 = cam1FromCam0.Inverse();

        int w0, h0;
        invWarp0.GetUVBounds(w0, h0);

        int w1, h1;
        invWarp1.GetUVBounds(w1, h1);
        Image<S> depth(w0, h0, 1);

        for (int i = 0; i < disparity.h; ++i)
        {
            for (int j = 0; j < disparity.w; ++j)
            {
                S d = disparity(i, j);
                
                S u0, v0;
                invWarp0.GetUV(j, i, u0, v0);

                if (u0 < S(0) || u0 + 1 >= w0 || v0 < S(0) || v0 + 1 >= h0)
                {
                    continue;
                }

                S u1, v1;
                invWarp1.InterpolateUV(j + d, i, u1, v1);

                if (u1 < S(0) || u1 + 1 >= w1 || v1 < S(0) || v1 + 1 >= h1)
                {
                    continue;
                }

                // Calibrations are for normalized coordinates
                u0 /= S(w0);
                v0 /= S(h0);
                u1 /= S(w1);
                v1 /= S(h1);

                Vec3<S> ray0;
                pinhole0.Unproject(u0, v0, ray0(0), ray0(1));
                ray0(2) = S(1);

                Vec3<S> ray1;
                pinhole1.Unproject(u1, v1, ray1(0), ray1(1));
                ray1(2) = S(1);

                Vec3<S> point;
                if (FindRayIntersection<S>(ray0, ray1, cam0FromCam1.t, point))
                {
                    // Store the z-value
                    depth(int(round(u0)), int(round(v0))) = point(2);
                }
            }
        }

        return depth;
    }
}