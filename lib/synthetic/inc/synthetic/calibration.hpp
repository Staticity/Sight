#pragma once

#include <linear/vec.hpp>
#include <calibration/cameramodel.hpp>

namespace sight
{
    /**
     * @brief Generates a list of 3D points on an XY grid.
     * 
     * @tparam S The floating point type
     * @param rows Number of grid rows
     * @param cols Number of grid columns
     * @param dx Distance in x between adjacent grid points
     * @param dy Distance in y between adjacent grid points
     * @return std::vector<Vec3<S>>
     */
    template <typename S>
    std::vector<Vec3<S>> GeneratePointGrid(
        const int rows,
        const int cols,
        const S dx,
        const S dy)
    {
        std::vector<Vec3<S>> points
        points.reserve(rows * cols);
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                points.push_back({i * dx, j * dy, S(0)});
            }
        }
        return points;
    }

    /**
     * @brief Projects 3D points in the camera into pixels in its image.
     * 
     * @tparam S The floating point type
     * @param model The camera model to project points
     * @param pointsInCamera List of 3D points in the camera's coordinate frame
     * @param filterPixels True to remove pixels outside the image's bounds.
     * @param maxW Maximum x-coordinate for imaged pixels
     * @param maxH Maximum y-coordinate for imaged pixels
     * @return std::vector<Vec2<S>> 
     */
    template <typename S>
    std::vector<Vec2<S>> ProjectPoints(
        const ICameraModel<S>& model,
        const std::vector<Vec3<S>>& pointsInCamera,
        bool filterPixels = false,
        const S maxW = S(1),
        const S maxH = S(1))
    {
        std::vector<Vec2<S>> pixels;
        pixels.reserve(pointsInCamera.size());
        
        if (filterPixels)
        {
            for (const auto& pt : pointsInCamera)
            {
                Vec2<S> pixel;
                model.Project(pt, pixel);

                if (pixel[0] > S(0) &&
                    pixel[1] > S(0) &&
                    pixel[0] < maxW &&
                    pixel[1] < maxH)
                {
                    pixels.push_back(std::move(pixel));
                }
            }
        }
        else
        {
            for (const auto& pt : pointsInCamera)
            {
                Vec2<S> pixel;
                model.Project(pt, pixel);
                pixels.push_back(std::move(pixel));
            }
        }

        return pixels;
    }
}
