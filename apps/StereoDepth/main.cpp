#include <iostream>
#include <filesystem>
#include <chrono>

#include <calibration/equidistantmodel.hpp>
#include <calibration/radial4model.hpp>
#include <geometric/rectify.hpp>
#include <geometric/polar_rectify.hpp>
#include <utils/ply_helpers.hpp>
#include <utils/timer.hpp>

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
    std::cout << t.DurationAndReset() << " seconds [reading]" << std::endl;
    
    // ComputeStereoRectificationTransforms()

    // Hard-coding the assignment of calibration right now..
    using S = double;
    Camera<S> cam0, cam1;
    cam0.CopyCamera(
        EquidistantModel<S> (
            380.81042871360756 / im0.w,
            380.81194179427075 / im0.h,
            510.29465304840727 / im0.w,
            514.3304630538506 / im0.h,
            0.010171079892421483,
            -0.010816440029919381,
            0.005942781769412756,
            -0.001662284667857643));
    cam0.Rt = {
        {
            -0.9995378259923383, 0.02917807204183088, -0.008530798463872679,
            0.007526588843243184, -0.03435493139706542, -0.9993813532126198,
            -0.029453096117288798, -0.9989836729399656, 0.034119442089241274
        },
        {
            0.047094247958417004, -0.04788273017221637, -0.0697294754693238
        }
    };

    cam1.CopyCamera(
        EquidistantModel<S> (
            379.2869884263036 / im1.w,
            379.26583742214524 / im1.h,
            505.5666703237407 / im1.w,
            510.2840961765407 / im1.h,
            0.01371679169245271,
            -0.015567360615942622,
            0.00905043103315326,
            -0.002347858896562788));
    cam1.Rt = {
        {
            -0.9995240747493029, 0.02986739485347808, -0.007717688852024281,
            0.008095979457928231, 0.01256553460985914, -0.9998882749870535,
            -0.02976708103202316, -0.9994748851595197, -0.0128013601698453
        },
        {
            -0.05374086123613335, -0.04648588412432889, -0.07333210787623645
        }
    };

    // Undistort the images and retrieve the new, associated pinhole models
    Image<uint8_t> u0, u1;
    PinholeModel<S> K0, K1;
    const S radius0 = cam0.model->ComputeProjectiveRadius();
    UndistortImage(im0, *cam0.model, u0, K0, radius0, radius0);
    UndistortImage(im1, *cam1.model, u1, K1, radius0, radius0);

    // After distortion, rectify the images using Pollefeys' polar coordinate transformation
    Image<uint8_t> r0, r1;
    InversePolarWarp<S> iw0, iw1;
    PolarRectifyImages(u0, K0, u1, K1, cam1.Rt * cam0.Rt.Inverse(), r0, iw0, r1, iw1);
    Image<S> disparity = DisparityFromStereoRectified<uint8_t, S>(r0, r1);

    std::vector<Vec3<S>> points;
    Image<S> depth = DepthFromDisparity(disparity, K0, iw0, K1, iw1, cam1.Rt * cam0.Rt.Inverse(), &points);
    WritePlyPointCloud("points.ply", points);    

    S minV, maxV;
    MinMax(disparity, minV, maxV);
    for (int i = 0; i < disparity.h; ++i)
        for (int j = 0; j < disparity.w; ++j)
            disparity(i, j) = (disparity(i, j) - minV) / (maxV - minV) * S(255);

    MinMax(depth, minV, maxV);
    for (int i = 0; i < depth.h; ++i)
        for (int j = 0; j < depth.w; ++j)
            depth(i, j) = (depth(i, j) - minV) / (maxV - minV) * S(255);

    std::cout << t.DurationAndReset() << " seconds [stereo]" << std::endl;

    const auto mu0 = Resize<uint8_t>(u0, .5f, .5f);
    const auto mu1 = Resize<uint8_t>(u1, .5f, .5f);
    const auto im02 = Resize<uint8_t>(im0, .5f, .5f);

    cv::imshow("Undistorted 0", mu0.ToOpenCV());
    cv::imshow("Undistorted 1", mu1.ToOpenCV());
    cv::imshow("Rect 0", r0.ToOpenCV());
    cv::imshow("Rect 1", r1.ToOpenCV());
    cv::imshow("Disparity", disparity.ToOpenCV());
    cv::imwrite("Disparity.png", disparity.ToOpenCV());
    cv::imshow("Depth", depth.ToOpenCV());
    cv::imwrite("Depth.png", depth.ToOpenCV());
    cv::imshow("Original", im02.ToOpenCV());
    cv::imwrite("PolarRectified0.png", r0.ToOpenCV());
    cv::imwrite("PolarRectified1.png", r1.ToOpenCV());
    cv::waitKey(0);
}