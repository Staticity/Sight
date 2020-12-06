#pragma once

#include <calibration/cameramodel.hpp>
#include <io/video_sources.hpp>
#include <utils/string_helpers.hpp>

#include <yaml-cpp/yaml.h>

#include <combaseapi.h>

#include <ctime>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

namespace sight
{
    template <typename S>
    struct CameraCalibrationIO
    {
        std::string identifier;
        std::time_t date;
        std::string dateStr;
        CameraModel<S> model;
    };

    struct CameraInfoIO
    {
        std::string identifier;
        std::string guid;
    };

    /*
    *  CalibrationStore contains the 'database' of camera calibrations
    *  across all devices seen on this (Windows) machine.
    *
    *  Internal directory structure:
    *    - <store directory>/
    *        - config.yml
    *        - <cameras dir>/
    *            - <camera id>/
    *                - calibration.yml
    *
    *
    *  'config.yml' format:
    *      cameras-dir: <insert>
    *
    *  'calibration.yml' format:
    *      identifier: <insert>
    *      date: <insert>
    *      model: <insert>
    *
    */
    class CalibrationStore
    {
    public:
        CalibrationStore() = default;
        CalibrationStore(const std::filesystem::path& storePath)
        {
            Load(storePath);
        }

        const std::filesystem::path& Path() const
        {
            return m_storePath;
        }

        bool Load(const std::filesystem::path& storePath)
        {
            // The path must be a directory
            if (!std::filesystem::is_directory(storePath))
            {
                // Attempt to create it before failing
                if (!std::filesystem::create_directories(storePath))
                {
                    return false;
                }
            }

            // Make sure there's a config file
            YAML::Node config;
            const auto configPath = storePath / "config.yml";
            if (!std::filesystem::exists(configPath))
            {
                // We'll create a default config file if it doesn't exist
                config["cameras-dir"] = "cameras";
                config["cameras"] = std::vector<CameraInfoIO>();
                
                std::ofstream out(configPath.string());
                out << config;

                return false;
            }
            else
            {
                config = YAML::LoadFile(configPath.string());
            }

            // Check if any fields are missing
            const std::vector<std::string> expectedFields = { "cameras-dir", "cameras" };
            if (std::any_of(
                expectedFields.begin(),
                expectedFields.end(),
                [&](auto& f) { return !config[f].IsDefined(); }))
            {
                return false;
            }

            // Store the fields
            m_storePath = storePath;
            m_config = config;

            m_cameraInfos = m_config["cameras"].as<std::vector<CameraInfoIO>>();
            for (int i = 0; i < m_cameraInfos.size(); ++i)
            {
                const auto& info = m_cameraInfos[i];
                m_cameraIdToIndex[info.identifier] = i;
            }

            return true;
        }

        template <typename S>
        CameraCalibrationIO<S> LoadCamera(int connectedCameraIdx) const
        {
            const auto source = GetVideoSource(connectedCameraIdx);
            if (!source.valid)
            {
                return {};
            }

            return LoadCamera<S>(source.identifier);
        }

        template <typename S>
        CameraCalibrationIO<S> LoadCamera(const std::string& identifier) const
        {
            const auto path = GetCameraCalibrationDirectory(identifier) / "calibration.yml";
            if (!std::filesystem::exists(path))
            {
                return {};
            }

            return YAML::LoadFile(path.string()).as<CameraCalibrationIO<S>>();
        }

        template <typename S>
        bool WriteCamera(int connectedCameraIdx, const CameraModel<S>& model)
        {
            const auto source = GetVideoSource(connectedCameraIdx);
            if (!source.valid)
            {
                return {};
            }
            return WriteCamera(source.identifier, model);
        }

        template <typename S>
        bool WriteCamera(const std::string& identifier, const CameraModel<S>& model)
        {
            const auto now = std::chrono::system_clock::now();
            const auto time = std::chrono::system_clock::to_time_t(now);
            const auto dateTm = std::gmtime(&time);

            CameraCalibrationIO<S> cal;
            cal.date = time;
            cal.dateStr = std::asctime(dateTm);
            cal.identifier = identifier;
            cal.model = model;

            YAML::Node node;
            node = cal;

            RegisterCamera(identifier);
            const auto dir = GetCameraCalibrationDirectory(identifier);
            const auto path = dir / "calibration.yml";

            if (!std::filesystem::exists(dir))
            {
                if (!std::filesystem::create_directories(dir))
                {
                    return false;
                }
            }

            std::ofstream calOut(path);
            std::ofstream configOut(m_storePath / "config.yml");
            calOut << node;
            configOut << m_config;

            return true;
        }


    private:

        void RegisterCamera(const std::string& identifier)
        {
            // Don't re-register
            if (m_cameraIdToIndex.count(identifier))
            {
                return;
            }

            GUID g;
            CoCreateGuid(&g);

            wchar_t guidWStr[256];
            StringFromGUID2(g, guidWStr, 256);

            const std::string guid = string_cast(guidWStr);
            CameraInfoIO info;
            info.identifier = identifier;
            info.guid = guid;

            m_config["cameras"].push_back(info);
            m_cameraInfos.push_back(info);
            m_cameraIdToIndex[info.identifier] = static_cast<int>(m_cameraInfos.size()) - 1;
        }

        std::optional<CameraInfoIO> FindCameraInfo(const std::string& identifier) const
        {
            if (!m_cameraIdToIndex.count(identifier))
            {
                return {};
            }
            return m_cameraInfos[m_cameraIdToIndex.at(identifier)];
        }

        std::string FindCameraGuid(const std::string& identifier) const
        {
            return FindCameraInfo(identifier).value_or(CameraInfoIO()).guid;
        }

        std::filesystem::path GetCameraCalibrationDirectory(const std::string& identifier) const
        {
            const auto camerasDir = m_config["cameras-dir"].as<std::string>();
            const auto deviceDir = FindCameraGuid(identifier);

            return m_storePath / camerasDir / deviceDir;
        }

        std::filesystem::path m_storePath;
        YAML::Node m_config;

        std::vector<CameraInfoIO> m_cameraInfos;
        std::map<std::string, int> m_cameraIdToIndex;
    };
}

namespace YAML
{
    template <>
    struct convert<sight::CameraInfoIO>
    {
        static Node encode(const sight::CameraInfoIO& cal)
        {
            Node node;
            node["identifier"] = cal.identifier;
            node["guid"] = cal.guid;
            return node;
        }

        static bool decode(const Node& node, sight::CameraInfoIO& cal)
        {
            const std::vector<std::string> fields = { "identifier", "guid" };
            for (const auto& field : fields)
            {
                if (!node[field].IsDefined())
                {
                    return false;
                }
            }

            cal.identifier = node["identifier"].as<std::string>();
            cal.guid = node["guid"].as<std::string>();
            return true;
        }
    };

    template <typename S>
    struct convert<sight::CameraCalibrationIO<S>>
    {
        static Node encode(const sight::CameraCalibrationIO<S>& cal)
        {
            Node node;
            node["identifier"] = cal.identifier;
            node["date"] = cal.date;
            node["date-friendly"] = cal.dateStr;
            node["model"] = cal.model;
            return node;
        }

        static bool decode(const Node& node, sight::CameraCalibrationIO<S>& cal)
        {
            const std::vector<std::string> fields = { "identifier", "date", "date-friendly", "model" };
            const auto isMissing = [&](auto& f) { return !node.IsDefined(); };
            if (std::any_of(
                    fields.begin(),
                    fields.end(),
                    isMissing))
            {
                return false;
            }

            cal.identifier = node[fields[0]].as<std::string>();
            cal.date = node[fields[1]].as<std::time_t>();
            cal.dateStr = node[fields[2]].as<std::string>();
            cal.model = node[fields[3]].as<sight::CameraModel<S>>();
            return true;
        }
    };
}