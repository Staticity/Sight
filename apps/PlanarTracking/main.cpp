#include <iostream>
#include <filesystem>

#include <Eigen/Eigen>
#include <optimization/single_camera.hpp>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    using namespace sight;
    using S = double;

    auto cam = PinholeModel<S>(S(.5), S(.5), S(.5), S(.5));
    SE3<S> pose = SE3<S>::Identity();
    OptimizeCameraIteration(std::vector<Frame<S>>(), cam, pose);
}
