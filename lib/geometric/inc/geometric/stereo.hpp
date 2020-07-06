#pragma once

#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <image/image_gen.hpp>
#include <calibration/pinholemodel.hpp>
#include <calibration/device.hpp>

namespace sight
{
    template <typename T, typename S>
    bool UndistortImage(
        const Image<T>& im,
        const CameraModel<S>* cam,
        Image<T>& outIm,
        PinholeModel<S>* outCam,
        S xRadius = std::numeric_limits<S>::max(),
        S yRadius = std::numeric_limits<S>::max())
    {
        if (outIm.w == 0 || outIm.h == 0)
        {
            outIm = Image<T>(im.w, im.h, im.c);
        }

        // Let's unproject all of the boundary pixels to
        // find the bounding rectangle of all the coordinates.
        const S invW = S(1);// / im.w;
        const S invH = S(1);// / im.h;
        S minX = std::numeric_limits<S>::max();
        S minY = std::numeric_limits<S>::max();
        S maxX = std::numeric_limits<S>::min();
        S maxY = std::numeric_limits<S>::min();
#if 0
        for (int i = 0; i < im.h; ++i)
        {
            // Question: Do I need to center these pixels..?
            //
            // Depends on how I want to do the camera calibration.
            // I assume that I would need to center it, so I will.
            S xz, yz;

            // Left column
            cam->Unproject(S(0.5) * invW, (S(i) + S(0.5)) * invH, xz, yz);
            minX = std::min(xz, minX);
            minY = std::min(yz, minY);
            maxX = std::max(xz, maxX);
            maxY = std::max(yz, maxY);

            // Right column
            cam->Unproject((S(im.w - 1) + S(0.5)) * invW, (S(i) + S(0.5)) * invH, xz, yz);
            minX = std::min(xz, minX);
            minY = std::min(yz, minY);
            maxX = std::max(xz, maxX);
            maxY = std::max(yz, maxY);
        }

        for (int j = 0; j < im.w; ++j)
        {
            S xz, yz;

            // Top row
            cam->Unproject((S(j) + S(0.5)) * invW, S(0.5) * invH, xz, yz);
            minX = std::min(xz, minX);
            minY = std::min(yz, minY);
            maxX = std::max(xz, maxX);
            maxY = std::max(yz, maxY);

            // Bottom row
            cam->Unproject((S(j) + S(0.5)) * invW, (S(im.h - 1) + S(0.5)) * invH, xz, yz);
            minX = std::min(xz, minX);
            minY = std::min(yz, minY);
            maxX = std::max(xz, maxX);
            maxY = std::max(yz, maxY);
        }
#else
        // We choose to unproject the entirety of the image
        // insetad of just the boundary, since some camera
        // models may not be able to unproject the boundaries
        // due to extreme distortion.
        //
        // We will slightly optimize this process by scanning
        // each scanline left-to-right until success, and vice-versa.
        for (int i = 0; i < im.h; ++i)
        {
            S xz, yz;

            int l;
            for (l = 0; l < im.w; ++l)
            {
                if (cam->Unproject((S(l) + S(0.5)) * invW, (S(i) + S(0.5)) * invH, xz, yz))
                {
                    minX = std::min(xz, minX);
                    minY = std::min(yz, minY);
                    maxX = std::max(xz, maxX);
                    maxY = std::max(yz, maxY);
                    break;
                }
            }

            for (int r = im.w - 1; r > l; --r)
            {
                if (cam->Unproject((S(r) + S(0.5)) * invW, (S(i) + S(0.5)) * invH, xz, yz))
                {
                    minX = std::min(xz, minX);
                    minY = std::min(yz, minY);
                    maxX = std::max(xz, maxX);
                    maxY = std::max(yz, maxY);
                    break;
                }
            }
        }
#endif

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

        const S fx = outIm.w / dx;
        const S fy = outIm.h / dy;
        const S cx = S(0.5) * invW - fx * minX;
        const S cy = S(0.5) * invH - fy * minY;

        outCam->Param(PinholeModel<S>::FX) = fx;
        outCam->Param(PinholeModel<S>::FY) = fy;
        outCam->Param(PinholeModel<S>::CX) = cx;
        outCam->Param(PinholeModel<S>::CY) = cy;
        outCam->Param(PinholeModel<S>::SKEW) = S(0);
        
        // Now, we must inverse sample the pixels in the original
        // image into the new image.
        auto view = PadView(im, 1);
        const S destInvW = S(1);// / outIm.w;
        const S destInvH = S(1);// / outIm.h;
        for (int i = 0; i < outIm.h; ++i)
        {
            T* row = outIm.row(i);
            for (int j = 0, k = 0; j < outIm.w; ++j)
            {
                bool ok = true;
                // Map back into metric space
                S xz, yz;
                ok &= outCam->Unproject(
                    (S(j) + S(.5)) * invW,
                    (S(i) + S(.5)) * invH,
                    xz,
                    yz);

                // Project into the original image
                S u, v;
                ok &= cam->Project(xz, yz, S(1), u, v, nullptr, nullptr);
                // u = (u * im.w) - S(.5);
                // v = (v * im.h) - S(.5);
                u = (u ) - S(.5);
                v = (v ) - S(.5);

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

        return true;
    }


    template <typename T, typename S>
    bool StereoRectify(
        const Image<T>& im0,
        const Camera<S>& cam0,
        const Image<T>& im1,
        const Camera<S>& cam1,
        Image<T>& rectIm0,
        Camera<S>& recCam0,
        Image<T>& rectIm1,
        Camera<S>& recCam1)
    {
        // This should create a pair of images, where
        // rectIm0.row(i) corresponds to rectIm1.row(i).
        //
        
        // Unproject all of the points in both images to z=1.
        //
        // Find the bounding rectangle of both sets of unprojected
        // points (metric rays).

        // When you rotate a camera, you rotate the plane by the opposite
        // transformation.
        //
        // Comptute half of that transformation to distribute to both
        // cameras.
        const SO3<S> R1FromR0 = cam1.Rt.R * cam0.Rt.R.Inverse();
        const SO3<S> R1FromR0 = SO3<S>::Exp(S(.5) * cam1FromCam0.Log());
        const auto& updateR0 = halfCam1FromCam0;
        const auto updateR1 = updateT0.Inverse();
        
        return false;
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
        assert(cam0.model.Name() == "pinhole");
        assert(cam1.model.Name() == "pinhole");

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