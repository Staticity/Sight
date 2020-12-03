#include <iostream>
#include <fstream>
#include <vector>

#include <calibration/radialtanmodel.hpp>
#include <geometric/rectify.hpp>
#include <image/image.hpp>
#include <io/video_device.hpp>

#include <Eigen/Eigen>

#include <opencv2/opencv.hpp>
#include <opencv2/core/utils/logger.hpp>

#include <yaml-cpp/yaml.h>
#include <io/cameramodel_io.hpp>

int main(int argc, char** argv)
{
	using namespace sight;
	using namespace std;

	//cout << cv::getBuildInformation() << endl;
	//cv::utils::logging::setLogLevel(cv::utils::logging::LOG_LEVEL_VERBOSE);

	sight::VideoDevice device(0);// R"(\\\\?\\usb#vid_32e4&pid_0144&mi_00#6&1e3f9547&1&0000#{65e8773d-8f56-11d0-a3b9-00a0c9223196}\\global)");

	auto node = YAML::LoadFile(R"(C:\Users\Jaime\source\repos\Sight\bin\elp_camera1.yml)");
	using S = double;
	auto model = node.as<RadialTanModel<S>>();

	//cv::VideoCapture cap(1, cv::CAP_MSMF);
	if (!device.IsOpen())
	{
		cerr << "Couldn't open the camera" << endl;
		return -1;
	}

	int width = static_cast<int>(device.GetProperty(cv::CAP_PROP_FRAME_WIDTH));
	int height = static_cast<int>(device.GetProperty(cv::CAP_PROP_FRAME_HEIGHT));
	auto start = device.GetProperty(cv::CAP_PROP_POS_MSEC);

	PinholeModel<double> pinhole = FindOptimalLinearCalibration(
		model,
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
