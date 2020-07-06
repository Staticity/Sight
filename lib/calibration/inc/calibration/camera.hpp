#pragma once

#include <lie/se3.hpp>
#include <calibration/cameramodel.hpp>

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

        SE3<S> Rt;
        CameraModel<S> model;
    };
}