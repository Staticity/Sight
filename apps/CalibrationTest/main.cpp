#include <iostream>
#include <fstream>
#include <vector>

#include <calibration/radialtanmodel.hpp>
#include <geometric/rectify.hpp>
#include <image/image.hpp>
#include <io/cameramodel_io.hpp>
#include <io/calibrationstore.hpp>
#include <io/video_device.hpp>
#include <optimization/camera.hpp>

#include <Eigen/Eigen>

#include <opencv2/opencv.hpp>
#include <opencv2/core/utils/logger.hpp>

#include <yaml-cpp/yaml.h>


int main(int argc, char** argv)
{
    using namespace sight;
    using namespace std;

    CalibrationStore calStore;
    calStore.Load("calibration_store");

    sight::VideoDevice device(0);
    if (!device.IsOpen())
    {
        cerr << "Couldn't open the camera" << endl;
        return -1;
    }

    auto cal = calStore.LoadCamera<double>(device.GetSourceInfo().identifier);
    auto& model = cal.model;

    using S = float;
    std::vector<OptimizeCameraIndex> indices;
    std::vector<Vec2<S>*> points2D;
    std::vector<Vec3<S>*> points3D;
    std::vector<SE3<S>*> devicesFromWorld;
    std::vector<SE3<S>*> camsFromDevice;
    std::vector<ICameraModel<S>*> cameras;
    sight::OptimizeCameraIteration<S>(indices, points2D, points3D, devicesFromWorld, camsFromDevice, cameras);

    int width = static_cast<int>(device.GetProperty(cv::CAP_PROP_FRAME_WIDTH));
    int height = static_cast<int>(device.GetProperty(cv::CAP_PROP_FRAME_HEIGHT));
    const auto start = device.GetProperty(cv::CAP_PROP_POS_MSEC);

    PinholeModel<double> pinhole = FindOptimalLinearCalibration(
        *model.Impl(),
        width, height,
        width, height);

    //cap.set(cv::CAP_PROP_FOURCC, cv::VideoWriter::fourcc('M', 'J', 'P', 'G'));

    double fps = device.GetProperty(cv::CAP_PROP_FPS);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    double exposure = 0.0;
    while (device.IsOpen())
    {
        std::cout << (device.GetProperty(cv::CAP_PROP_POS_MSEC) - start) / 1e3 << std::endl;
        Image<uint8_t> frame = device.ReadFrame();

        //Image<uint8_t> image = Image<uint8_t>::FromOpenCV(frame);

        //Image<uint8_t> result;
        //UndistortImage(image, camera, pinhole, result);

        cv::imshow("Original", Image<uint8_t>::ToOpenCV(frame));
        //cv::imshow("Undistorted", Image<uint8_t>::ToOpenCV(result));

        char key;
        if (((key = cv::waitKey(1)) & 0xff) == 27) break;
    }

}
