#pragma once

#include <calibration/camera.hpp>

namespace sight
{
    // Contains a collection of relevant sensor information
    // for performing computer vision tasks -- such as their
    // calibrations in 'device space' containing poses and
    // other modeling/correction terms.
    template <typename S>
    struct Device
    {
        std::vector<Camera<S>> cameras;
    };
}