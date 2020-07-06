#include <iostream>
#include <filesystem>
#include <chrono>
#include <utils/timer.hpp>

#include <geometric/stereo.hpp>
#include <calibration/equidistantmodel.hpp>
#include <calibration/radial4model.hpp>

#include <opencv2/opencv.hpp>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    using namespace sight;

    const fs::path exe = __FILE__;
    const std::string folder = (argc > 1) ? argv[1] : "mav";
    const std::string im0Path = (exe.parent_path() / "inputs" / folder / "image0.png").string();
    const std::string im1Path = (exe.parent_path() / "inputs" / folder / "image1.png").string();

    Timer t;
    const auto im0 = ToGrayscale(Image<uint8_t>::Read(im0Path));
    const auto im1 = ToGrayscale(Image<uint8_t>::Read(im1Path));
    std::cout << t.DurationAndReset() << " seconds reading" << std::endl;
    
    using S = double;
    const EquidistantModel<S> cam0(
        380.81042871360756,
        380.81194179427075,
        510.29465304840727,
        514.3304630538506,
        0.010171079892421483,
        -0.010816440029919381,
        0.005942781769412756,
        -0.001662284667857643);

    Image<uint8_t> undistorted;
    PinholeModel<S> K;
    const bool success = UndistortImage(im0, &cam0, undistorted, &K, 3.0, 3.0);
    std::cout << t.DurationAndReset() << " seconds reading" << std::endl;

    const auto undistorted2 = Resize<uint8_t>(undistorted, .5f, .5f);
    const auto im02 = Resize<uint8_t>(im0, .5f, .5f);

    cv::imshow("Undistorted", undistorted2.ToOpenCV());
    cv::imshow("Original", im02.ToOpenCV());
    cv::waitKey(0);
}