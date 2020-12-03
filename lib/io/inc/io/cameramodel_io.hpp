#pragma once

#include <calibration/cameramodelfactory.hpp>
#include <yaml-cpp/node/convert.h>

#define GENERATE_CONVERSION(Model) \
    template <typename S> \
    struct convert<Model>\
    {\
        static Node encode(const sight::ICameraModel<S>& model)\
        {\
            return convert<sight::ICameraModel<S>>::encode(model);\
        }\
        static bool decode(const Node& node, sight::ICameraModel<S>& model)\
        {\
        return convert<sight::ICameraModel<S>>::decode(node, model);\
        }\
    };\

namespace YAML
{
    template <typename S>
    struct convert<sight::CameraModel<S>>
    {
        static bool ValidateFields(const Node& node)
        {
            std::vector<std::string> fields = { "name", "params" };
            for (const auto& field : fields)
            {
                if (!node[field].IsDefined())
                {
                    return false;
                }
            }
            return true;
        }

        static Node encode(const sight::CameraModel<S>& model)
        {
            std::vector<S> params(model.NumParams());
            for (int i = 0; i < params.size(); ++i)
            {
                params[i] = model.Param(i);
            }

            Node node;
            node["name"] = model.Name();
            node["params"] = params;
            return node;
        }

        static bool decode(const Node& node, sight::CameraModel<S>& model)
        {
            if (!ValidateFields(node))
            {
                return false;
            }

            const auto name = node["name"].as<std::string>();
            const auto params = node["params"].as<std::vector<S>>();

            // Ensure that our YAML file is up to date the
            // camera model factory.
            model = sight::CreateCameraModel<S>(name);
            if (params.size() != model.NumParams())
            {
                return false;
            }

            // Load each parameter into the model
            for (int i = 0; i < model.NumParams(); ++i)
            {
                model.Param(i) = params[i];
            }

            return true;
        }
    };

    template <typename S>
    struct convert<sight::ICameraModel<S>>
    {
        static Node encode(const sight::ICameraModel<S>& model)
        {
            return convert<sight::CameraModel<S>>::encode(model->Clone());
        }

        static bool decode(const Node& node, sight::ICameraModel<S>& model)
        {
            return convert<sight::CameraModel<S>>::decode(node, model->Clone());
        }

    };

}
