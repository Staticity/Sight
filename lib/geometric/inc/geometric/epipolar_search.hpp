#pragma once

#include <calibration/camera>
#include <linear/vec.hpp>

namespace sight
{
    
    template <typename T>
    bool EpipolarSearch(
        const Camera<S>& cam0,
        const Camera<S>& cam1,
        const Vec2<S>& uv,
        S& depth,
        S minDepth = std::numeric_limits<S>::epsilon(),
        S maxDepth = std::numeric_limits<S>::max())
    {
        const SE3<S> cam1FromCam0 = cam1.Rt * cam0.Rt.Inverse();

        Vec3<S> ray0;
        if (!cam0.model->Unproject(uv, ray0))
        {
            return false;
        }

        // As I change z, how will the uv change?
        
        S Jxyz[6];
        Vec2<S> uv_unused;
        cam0.model->Project(ray0, uv_unused, nullptr, Jxyz);
    }
}