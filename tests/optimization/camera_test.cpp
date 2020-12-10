#include <catch2/catch.hpp>

#include <calibration/cameramodel.hpp>
#include <io/cameramodel_io.hpp>
#include <optimization/camera.hpp>

#include<yaml-cpp/yaml.h>

#include <iostream>
#include <filesystem>

TEST_CASE("Optimize single camera, single frame")
{
    using namespace sight;
    using S = double;

    CameraModel<S> camera(new PinholeModel<S>(.5, .5, .5, .5));
    const SE3<S> deviceFromWorld = SE3<S>::Exp({S(.1), S(-.35), S(1.22), S(0.01), S(2.2), S(3.3)});
    const SE3<S> camFromDevice = SE3<S>::Identity();

    std::vector<Vec3<S>> xyzs;
    std::vector<Vec2<S>> uvs;
    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            xyzs.push_back({ S(i), S(j), S(1) });
            uvs.emplace_back();
            camera.Project(camFromDevice * deviceFromWorld * xyzs.back(), uvs.back());
        }
    }

    std::vector<OptimizeCameraIndex> indices;
    std::vector<Vec3<S>*> p_xyzs;
    std::vector<Vec2<S>*> p_uvs;
    for (int i = 0; i < xyzs.size(); ++i)
    {
        OptimizeCameraIndex index;
        index.index2D = i;
        index.index3D = i;
        index.indexCam = 0;
        index.indexSE3 = 0;

        indices.push_back(index);
        p_xyzs.push_back(&xyzs[i]);
        p_uvs.push_back(&uvs[i]);
    }

    const SE3<S> deviceNoise = SE3<S>::Exp({ S(.005), S(-.012), S(.001), S(-.01), S(.02), S(-.03) });

    CameraModel<S> camera_test(camera.Clone());
    SE3<S> deviceFromWorld_test = deviceNoise * deviceFromWorld;
    SE3<S> camFromDevice_test = camFromDevice;
    
    const int numIterations = 100;
    bool optimizationConverged = false;
    for (int i = 0; i < numIterations; ++i)
    {
        std::vector<SE3<S>*> devicesFromWorld = { &deviceFromWorld_test };
        std::vector<SE3<S>*> camsFromDevice = { &camFromDevice_test };
        std::vector<ICameraModel<S>*> cameras = { camera_test.Impl() };

        const S res = OptimizeCameraIteration<S>(
            indices,
            p_uvs,
            p_xyzs,
            devicesFromWorld,
            camsFromDevice,
            cameras,
            S(1.0));

        if (res < std::numeric_limits<S>::epsilon())
        {
            optimizationConverged = true;
            break;
        }
    }

    // Only requirement is that the residual goes
    // down to 0. It's possible there are gauge
    // freedoms between pose and model, so we
    // cannot compare those.
    REQUIRE(optimizationConverged);
}


TEST_CASE("Optimize dual camera, multiple frames")
{
    using namespace sight;
    using S = double;

    CameraModel<S> camera(new PinholeModel<S>(.5, .5, .5, .5));
    CameraModel<S> camera2(camera);
    const SE3<S> deviceFromWorld0 = SE3<S>::Exp({ S(.1), S(-.35), S(1.22), S(0.01), S(2.2), S(3.3) });
    const Vec<S, 6> deviceFromWorldVelocity = { S(.010), S(-.010), S(.001), S(.01), S(.02), S(.03) };
    SE3<S> camFromDevice = SE3<S>::Identity();
    SE3<S> camFromDevice2 = SE3<S>::Exp({ S(.001), S(.005), S(0.0001), S(-.01), S(.00015), S(-.0023) });

    std::vector<OptimizeCameraIndex> indices;
    std::vector<Vec2<S>*> p_uvs;
    std::vector<Vec3<S>*> p_xyzs;
    std::vector<SE3<S>*> devicesFromWorld;
    std::vector<SE3<S>*> camsFromDevice = { &camFromDevice, &camFromDevice2 };
    std::vector<ICameraModel<S>*> cameras = { camera.Impl(), camera2.Impl() };

    std::vector<Vec3<S>> xyzs(100);
    std::vector<Vec2<S>> uvs(100);

    for (int ti = 0; ti < 10; ++ti)
    {
        const double t = ti * .1;

        const auto deviceFromWorld = SE3<S>::Exp(deviceFromWorldVelocity * t) * deviceFromWorld0;
        devicesFromWorld.push_back(new SE3<S>(deviceFromWorld));

        for (int i = 0, k = 0; i < 10; ++i)
        {
            for (int j = 0; j < 10; ++j, ++k)
            {
                auto& index = indices.emplace_back();
                index.index2D = int(p_uvs.size());
                index.index3D = int(p_xyzs.size());
                index.indexCam = 0;
                index.indexSE3 = ti;

                xyzs[k] = { S(i), S(j), S(1) };
                camera.Project(camFromDevice * deviceFromWorld * xyzs[k], uvs[k]);
            }
        }
    }

    for (int i = 0; i < xyzs.size(); ++i)
    {
        p_xyzs.push_back(&xyzs[i]);
        p_uvs.push_back(&uvs[i]);
    }

    // Add noise to parameters
    for (int i = 0; i < camera.NumParams(); ++i)
    {
        camera.Param(i) += S(1e-2);
    }

    // Add noise to poses
    for (int i = 0; i < camsFromDevice.size(); ++i)
    {
        auto& Rt = *camsFromDevice[i];
        const S mult = (i % 2 == 0) ? S(.001) : S(-.001);
        Rt = SE3<S>::Exp(Vec<S, 6>(i * mult)) * Rt;
    }

    for (int i = 0; i < devicesFromWorld.size(); ++i)
    {
        auto& Rt = *devicesFromWorld[i];
        Vec<S, 6> v;
        for (int j = 0; j < 6; ++j)
        {
            v[j] = S(.001) * (j % 2) ? S(1) : S(-1);
        }
        Rt = SE3<S>::Exp(v) * Rt;
    }

    const int numIterations = 100;
    bool optimizationConverged = false;
    for (int i = 0; i < numIterations; ++i)
    {
        const S res = OptimizeCameraIteration<S>(
            indices,
            p_uvs,
            p_xyzs,
            devicesFromWorld,
            camsFromDevice,
            cameras,
            S(1.0));

        if (res < std::numeric_limits<S>::epsilon())
        {
            optimizationConverged = true;
            break;
        }
    }

    for (auto& deviceFromWorld : devicesFromWorld)
    {
        delete deviceFromWorld;
    }

    // Only requirement is that the residual goes
    // down to 0. It's possible there are gauge
    // freedoms between pose and model, so we
    // cannot compare those.
    REQUIRE(optimizationConverged);
}
