#include <iostream>
#include <filesystem>
#include <chrono>
#include <opencv2/opencv.hpp>
#include <utils/timer.hpp>

#include <features/optical_flow.hpp>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    using namespace sight;

    const fs::path exe = __FILE__;
    const std::string folder = (argc > 1) ? argv[1] : "warehouse";
    const std::string im0Path = (exe.parent_path() / "inputs" / folder / "image0.png").string();
    const std::string im1Path = (exe.parent_path() / "inputs" / folder / "image1.png").string();

    const auto im0 = ToGrayscale(Image<uint8_t>::Read(im0Path));
    const auto im1 = ToGrayscale(Image<uint8_t>::Read(im1Path));

    // {
    //     auto mat = Sobel3x3_Horizontal<float>(PadView(im0, 1)).ToOpenCV();
    //     cv::imshow("SobelX", mat);
    // }
    // {
    //     auto mat = Sobel3x3_Vertical<float>(PadView(im0, 1)).ToOpenCV();
    //     cv::imshow("SobelY", mat);
    // }
    // cv::waitKey(0);

    std::vector<Flow<float>> flows;
    OpticalFlow_KLT(im0, im1, flows, 10, 20, 100, 10);
 
    auto mat = ToRGB(im0).ToOpenCV();
    for (const auto& flow : flows)
    {
        const cv::Scalar color = (flow.valid) ? cv::Scalar(0, 255, 0) : cv::Scalar(0, 0, 255);
        const cv::Point2d center = { flow.pos[0], flow.pos[1] };
        const cv::Point2d offset = { flow.velocity[0], flow.velocity[1] };
        // cv::circle(mat, center, flow.radius, color);
        cv::arrowedLine(mat, center, center + 3 * offset, color);
    }
    cv::imshow("Flow", mat);
    cv::waitKey(0);
}