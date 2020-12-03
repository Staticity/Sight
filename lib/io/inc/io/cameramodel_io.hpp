#pragma once

#include <calibration/cameramodelfactory.hpp>
#include <yaml-cpp/node/convert.h>

#define GENERATE_CONVERSION(Model) \
    template <typename S> \
    struct convert<Model>\
    {\
        static Node encode(const sight::CameraModel<S>& model)\
        {\
            return convert<sight::CameraModel<S>>::encode(model);\
        }\
        static bool decode(const Node& node, sight::CameraModel<S>& outModel)\
        {\
        return convert<sight::CameraModel<S>>::decode(node, outModel);\
        }\
    };\

namespace YAML
{
    template <typename S>
    struct convert<sight::CameraModel<S>>
    {
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

        static bool decode(const Node& node, sight::CameraModel<S>& outModel)
        {
            if (!ValidateFields(node))
            {
                return false;
            }

            const auto name = node["name"].as<std::string>();
            const auto params = node["params"].as<std::vector<S>>();

            // Ensure that our YAML file is up to date the
            // camera model factory.
            auto model = sight::CreateCameraModel(name);
            if (model->Name() != rhs.Name() ||
                params.size() != model->NumParams())
            {
                return false;
            }

            // Instead of loading directly into the output model,
            // load into the newly created model.
            for (int i = 0; i < model->NumParams(); ++i)
            {
                model->Param(i) = params[i];
            }

            outModel.LoadModel(*model);

            return true;
        }
    };

    GENERATE_CONVERSION(sight::EquidistantModel<S>);
    GENERATE_CONVERSION(sight::PinholeModel<S>);
    GENERATE_CONVERSION(sight::Radial4Model<S>);
    GENERATE_CONVERSION(sight::RadialTanModel<S>);
}
