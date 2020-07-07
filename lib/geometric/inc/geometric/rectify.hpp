#pragma once

#include <random>
#include <time.h>

#include <Eigen/Eigen>

#include <calibration/pinholemodel.hpp>
#include <calibration/device.hpp>
#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <image/image_gen.hpp>
#include <geometric/homography_ops.hpp>

namespace sight
{
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
    Image<S> DepthFromStereoRectified(
        const Image<T>& im0,
        const CameraModel<S>& cam0,
        const Image<T>& im1,
        const CameraModel<S>& cam1,
        const int windowRadius = 2,
        const float maxDisparityPerc = .05f,
        const bool refineSubpixel = true)
    {
        assert(im0.c == 1);
        assert(im1.c == 1);
        assert(cam0->model.Name() == "pinhole");
        assert(cam1->model.Name() == "pinhole");

        // Since the images are rectified, we only need to perform
        // searches along their horizontal epipolar lines.
        const int maxDisparity = int(round(im1.w * maxDisparity));
        const int disparityRadius = (maxDisparity + 1) / 2;

        // For each pixel, form a patch around it and search along the
        // epipolar line to find a patch which matches most.
        Image<S> depth, disparity;

        for (int y = 0; y < im0.h; ++i)
        {
            const T* row0 = im0.row(y);
            const T* row1 = im1.row(y);
            for (int x = 0; x < im0.w; ++x, row++)
            {
                
            }
        }

        return Image<S>();
    }

}