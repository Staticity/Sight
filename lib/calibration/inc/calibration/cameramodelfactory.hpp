#pragma once

#include <calibration/cameramodel.hpp>
#include <calibration/allcameramodels.hpp>
#include <memory>

namespace sight
{
    template <typename S>
    std::unique_ptr<CameraModel<S>> CreateCameraModel(const std::string s)
    {
        if (Radial4Model<S>::ModelName() == s)
        {
            return new Radial4Model<S>();
        }
        if (RadialTanModel<S>::ModelName() == s)
        {
            return new RadialTanModel<S>();
        }
        if (PinholeModel<S>::ModelName() == s)
        {
            return new PinholeModel<S>();
        }
        if (EquidistantModel<S>::ModelName() == s)
        {
            return new EquidistantModel<S>();
        }

        return nullptr;
    }
}
