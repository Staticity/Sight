#include <io/video_device.hpp>
#include <io/video_sources.hpp>

namespace sight
{
    VideoDevice::VideoDevice(int index)
    {
        auto sources = EnumerateVideoSources();
        if (index < sources.size())
        {
            m_sourceInfo = sources[index];
            m_capture.open(index, cv::CAP_MSMF);
        }
    }

    VideoDevice::VideoDevice(const std::string& devicePath)
    {
        auto sources = EnumerateVideoSources();

        int index = 0;
        while (index < sources.size())
        {
            if (sources[index].identifier == devicePath)
            {
                break;
            }
            ++index;
        }

        if (index < sources.size())
        {
            m_sourceInfo = sources[index];
            m_capture.open(index, cv::CAP_MSMF);
        }
    }

    bool VideoDevice::IsOpen() const
    {
        return m_capture.isOpened();
    }

    const VideoSourceInfo& VideoDevice::GetSourceInfo() const
    {
        return m_sourceInfo;
    }

    bool VideoDevice::SetProperty(int prop, double value)
    {
        return m_capture.set(prop, value);
    }

    double VideoDevice::GetProperty(int prop)
    {
        return m_capture.get(prop);
    }

    Image<uint8_t> VideoDevice::ReadFrame()
    {
        cv::Mat frame;
        m_capture >> frame;
        return Image<uint8_t>::FromOpenCV(frame);
    }
}