#pragma once

#include <vector>

#include <Eigen/Core>

#include <calibration/allcameramodels.hpp>
#include <calibration/camera.hpp>
#include <linear/vec.hpp>
#include <lie/se3.hpp>

namespace sight
{
    template <typename S>
    struct Correspondence
    {
        Vec3<S> point;
        Vec2<S> pixel;
    };

    template <typename S>
    struct Frame
    {
        std::vector<Correspondence<S>> correspondences;
    };

    template <typename S>
    Eigen::Matrix3<S> Skew(S x, S y, S z)
    {
        Eigen::Matrix3<S> skew;
        skew <<
            S(0),   -z,   y,
               z, S(0),  -x,
              -y,    x, S(0);
        return skew;
    }

    template <typename S>
    S OptimizeCameraIteration(
        const std::vector<Frame<S>>& frames,
        ICameraModel<S>& cam,
        SE3<S>& camFromFrame,
        S lambda = S(0))
    {
        using Eigen::Dynamic;
        using Eigen::RowMajor;
        using Eigen::Matrix;
        using Eigen::Vector;
        using Eigen::Map;

        S totalResidual = S(0);

        // First N parameters belong to the camera
        // Last 6 parameters belong to the pose
        Vector<S, Dynamic> g;
        Matrix<S, Dynamic, Dynamic, RowMajor> JTJ;

        const int N_cam = cam.NumParams();

        Matrix<S, 2, Dynamic, RowMajor> J;
        Matrix<S, Dynamic, Dynamic, RowMajor> J_uv_cam;
        Matrix<S, 2, 6> J_uv_pose;
        Matrix<S, 2, 3, RowMajor> J_uv_xyz;
        Matrix<S, 3, 6> J_xyz_pose;

        J_uv_cam.resize(2, N_cam);
        JTJ.resize(N_cam + 6, N_cam + 6);
        JTJ.setZero();
        g.resize(N_cam + 6);
        g.setZero();

        for (const auto& frame : frames)
        {
            for (const auto& pair : frame.correspondences)
            {
                const auto& pointInModel = pair.point;
                const auto& uvInCam = pair.pixel;

                const auto pointInCam = camFromFrame * pointInModel;

                Vec2<S> uvPredicted;
                cam.Project(pointInCam, uvPredicted, J_uv_cam.data(), J_uv_xyz.data());
                
                const auto residual = uvPredicted - uvInCam;
                totalResidual += residual.SquaredNorm();

                // Compute J_uv_pose
                J_xyz_pose.block<3, 3>(0, 0) = Skew(-pointInCam(0), -pointInCam(1), -pointInCam(2));
                J_xyz_pose.block<3, 3>(3, 0) = Eigen::Matrix3<S>::Identity();
                J_uv_pose = J_uv_xyz * J_xyz_pose;

                J.block(0, 0, 2, N_cam) = J_uv_cam;
                J.block<2, 6>(N_cam, 0) = J_uv_pose;

                const auto JT = J.transpose();
                JTJ += JT * J;

                g += JT * Vector<S, 2>(residual(0), residual(1));
            }
        }
        
        // Dampen the Hessian
        Matrix<S, Dynamic, Dynamic, RowMajor> H = JTJ;
        for (int i = 0; i < H.rows(); ++i)
        {
            H(i, i) += lambda;
        }

        // Solve
        const Vector<S, Dynamic> x = H.llt().solve(-g);

        // Update camera intrinsics
        for (int i = 0; i < N_cam; ++i)
        {
            cam.Param(i) += x(i);
        }

        // Update camera extrinsics
        Vec<S, 6> log;
        for (int i = 0; i < 6; ++i)
        {
            log[i] = x(N_cam + i);
        }

        camFromFrame = SE3<S>::Exp(log) * camFromFrame;

        return totalResidual;
    }
}
