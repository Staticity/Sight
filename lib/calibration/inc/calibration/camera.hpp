#pragma once

#include <lie/se3.hpp>
#include <calibration/cameramodel.hpp>

#include <memory>

namespace sight
{
    template <typename S>
    class Camera
    {
    public:
        Camera()
            : Rt(SE3<S>::Identity())
            , model()
        {}

        void CopyCamera(const CameraModel<S>& cam)
        {
            model.reset(cam.Clone());
        }

        SE3<S> Rt;
        std::shared_ptr<CameraModel<S>> model;
    };
}