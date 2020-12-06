#pragma once

#include <calibration/cameramodel.hpp>
#include <calibration/allcameramodels.hpp>
#include <memory>

namespace sight
{
    template <typename S>
    CameraModel<S> CreateCameraModel(const std::string s)
    {
        if (EquidistantModel<S>::ModelName() == s)
        {
            return new EquidistantModel<S>();
        }
        if (PinholeModel<S>::ModelName() == s)
        {
            return new PinholeModel<S>();
        }
        if (Radial4Model<S>::ModelName() == s)
        {
            return new Radial4Model<S>();
        }
        if (RadialTanModel<S>::ModelName() == s)
        {
            return new RadialTanModel<S>();
        }

        return nullptr;
    }
}
