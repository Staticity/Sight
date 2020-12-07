#pragma once

#include <calibration/allcameramodels.hpp>
#include <calibration/camera.hpp>
#include <linear/vec.hpp>
#include <lie/se3.hpp>
#include <utils/eigentypes.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include <numeric>
#include <vector>

namespace sight
{
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

    struct OptimizeCameraIndex
    {
        int index2D{ -1 };
        int index3D{ -1 };
        int indexSE3{ -1 };
        int indexCam{ -1 };
    };

    //   | cameras  poses
    // -----------------
    // c |
    //   |
    // p |
    //   |
    //
    // Every camera is independent of one
    // another. Every pose is independent
    // of another.

    template <typename S>
    S OptimizeCameraIteration(
        const std::vector<OptimizeCameraIndex>& indices,
        const std::vector<Vec2<S>*>& points2D,
        const std::vector<Vec3<S>*>& points3D,
        std::vector<SE3<S>*>& camFromWorlds,
        std::vector<ICameraModel<S>*> cameras,
        S lambda = S(0))
    {
        using Eigen::Dynamic;
        using Eigen::RowMajor;
        using Eigen::Matrix;
        using Eigen::Vector;
        using Eigen::Map;

        S totalSqResidual = S(0);

        int numCameraParams = 0;
        std::vector<int> cameraIndices;
        for (int i = 0; i < cameras.size(); ++i)
        {
            cameraIndices.push_back(numCameraParams);
            numCameraParams += cameras[i]->NumParams();
        }

        int numFrameParams = 0;
        std::vector<int> frameIndices;
        for (int i = 0; i < camFromWorlds.size(); ++i)
        {
            frameIndices.push_back(numCameraParams + numFrameParams);
            numFrameParams += 6;
        }

        const int numTotalParams = numCameraParams + numFrameParams;
        Matrix<S, Dynamic, Dynamic, RowMajor> JTJ;
        JTJ.resize(numTotalParams, numTotalParams);
        JTJ.setZero();

        Vector<S, Dynamic> g;
        g.resize(numTotalParams);
        g.setZero();

        // Convenience lambdas
        auto GetHessianBlock_CC = [&](int i) {
            const int n = cameras[i]->NumParams();
            return JTJ.block(cameraIndices[i], cameraIndices[i], n, n);
        };

        auto GetHessianBlock_CF = [&](int i, int j) {
            const int n = cameras[i]->NumParams();
            return JTJ.block(cameraIndices[i], frameIndices[j], n, 6);
        };

        auto GetHessianBlock_FC = [&](int i, int j) {
            const int n = cameras[j]->NumParams();
            return JTJ.block(frameIndices[i], cameraIndices[j], 6, n);
        };

        auto GetHessianBlock_FF = [&](int i) {
            return JTJ.block<6, 6>(frameIndices[i], frameIndices[i]);
        };

        auto GetGradientBlock_C = [&](int i) {
            return g.block(cameraIndices[i], 0, cameras[i]->NumParams(), 1);
        };

        auto GetGradientBlock_F = [&](int i) {
            return g.block<6, 1>(frameIndices[i], 0);
        };

        // Single storage for common jacobians
        EigenList<Matrix<S, Dynamic, Dynamic, RowMajor>> J_uv_cams(cameras.size());
        for (int i = 0; i < cameras.size(); ++i)
        {
            J_uv_cams[i].resize(2, cameras[i]->NumParams());
        }

        Vec2<S> uv;
        Matrix<S, 2, 6> J_uv_pose;
        Matrix<S, 2, 3, RowMajor> J_uv_xyz;
        Matrix<S, 3, 6> J_xyz_pose;

        for (const auto& index : indices)
        {
            const auto& point2D = *points2D[index.index2D];
            const auto& point3D = *points3D[index.index3D];
            const auto& camFromWorld = *camFromWorlds[index.indexSE3];
            const auto& camera = cameras[index.indexCam];

            auto& J_uv_cam = J_uv_cams[index.indexCam];

            const auto Rp = camFromWorld.R * point3D;
            const auto pointInCam = Rp + camFromWorld.t;

            camera->Project(pointInCam, uv, J_uv_cam.data(), J_uv_xyz.data());

            const auto residual = uv - point2D;
            totalSqResidual += residual.SquaredNorm();

            J_xyz_pose.block<3, 3>(0, 0) = Skew(-Rp(0), -Rp(1), -Rp(2));
            J_xyz_pose.block<3, 3>(3, 0) = Eigen::Matrix3<S>::Identity();
            J_uv_pose = J_uv_xyz * J_xyz_pose;

            const auto& camIdx = index.indexCam;
            const auto& se3Idx = index.indexSE3;
            GetHessianBlock_CC(camIdx) = J_uv_cam.transpose() * J_uv_cam;
            GetHessianBlock_FF(se3Idx) = J_uv_pose.transpose() * J_uv_pose;
            GetHessianBlock_CF(camIdx, se3Idx) = J_uv_cam.transpose() * J_uv_pose;
            GetHessianBlock_CF(se3Idx, camIdx) = GetHessianBlock_CF(camIdx, se3Idx);

            const Vector<S, 2> eResidual(residual(0), residual(1));
            GetGradientBlock_C(camIdx) = J_uv_cam.transpose() * eResidual;
            GetGradientBlock_F(se3Idx) = J_uv_pose.transpose() * eResidual;
        }

        // [ A  B] [x0] = [y0]
        // [B^T C] [x1]   [y1]
        //
        //  A  * x0 + B * x1 = y0
        // B^T * x0 + C * x1 = y1
        //
        // x0 = (A - B * C^-1 * B^T)^-1 * (y0 - B * C^-1 * y1)
        // x1 = C^-1 * (y1 - B^T * x0)
        const auto& A = JTJ.block(0, 0, numCameraParams, numCameraParams);
        const auto& B = JTJ.block(0, numCameraParams, numCameraParams, numFrameParams);
        const auto& C = JTJ.block(numCameraParams, numCameraParams, numFrameParams, numFrameParams);

        const auto& y0 = g.block(0, 0, numCameraParams, 1);
        const auto& y1 = g.block(numCameraParams, 0, numFrameParams, 1);

        // Since C is block diagonal, its inverse is
        // also block diagonal with the inverse of each block.
        Matrix<S, Dynamic, Dynamic> Cinv;
        Cinv.resize(numFrameParams, numFrameParams);
        Cinv.setZero();
        for (int i = 0; i < camFromWorlds.size(); ++i)
        {
            const int j = frameIndices[i];
            Cinv.block<6, 6>(j, j) = GetHessianBlock_FF(i).fullPivLu().inverse();
        }

        // Compute x0 using full inverse of C
        const auto M = A - B * Cinv * B.transpose();
        const auto b = y0 - B * Cinv * y1;
        Vector<S, Dynamic> x0 = M.ldlt().solve(b);

        // Solve for x1 by using LDLT solve per frame, given x0
        Vector<S, Dynamic> x1;
        x1.resize(numFrameParams);
        x1.setZero();
        const auto y1_prime = y1 - B.transpose() * x0;
        for (int i = 0; i < camFromWorlds.size(); ++i)
        {
            const int j = frameIndices[i];
            auto& x1_i = x1.block<6, 1>(j, 0);
            const auto& y1p_i = y1_prime.block<6, 1>(j, 0);

            x1_i = Cinv.block<6, 6>(j, j).ldlt().solve(y1p_i);
        }

        // Update all of the camera parameters
        for (int i = 0; i < cameras.size(); ++i)
        {
            auto& camera = cameras[i];
            const int index = cameraIndices[i];
            for (int j = 0; j < camera->NumParams(); ++j)
            {
                camera->Param(j) += x0(index + j);
            }
        }

        // Update all of the pose parameters
        for (int i = 0; i < camFromWorlds.size(); ++i)
        {
            auto& pose = *camFromWorlds[i];
            const int index = frameIndices[i];

            const Eigen::Vector3<S> dr = x1.block<3, 1>(index, 0);
            const Eigen::Vector3<S> dt = x1.block<3, 1>(index + 3, 0);

            const auto dR = SO3<S>::Exp({ dr(0), dr(1), dr(2) });
            pose.R = dR * pose.R;
            pose.t += { dt(0), dt(1), dt(2) };
        }

        return totalSqResidual;
    }
}
