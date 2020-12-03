#pragma once

#include <string>

#include <image/image.hpp>
#include <io/video_sources.hpp>
#include <opencv2/videoio/videoio.hpp>

namespace sight
{
    class VideoDevice
    {
    public:
        VideoDevice(int index);
        VideoDevice(const std::string& devicePath);

        bool IsOpen() const;
        const VideoSourceInfo& GetSourceInfo() const;

        bool SetProperty(int prop, double value);
        double GetProperty(int prop);

        Image<uint8_t> ReadFrame();

    private:
        VideoSourceInfo m_sourceInfo;
        cv::VideoCapture m_capture;
    };
}
