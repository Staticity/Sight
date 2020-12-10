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
    Eigen::Matrix3<S> Skew(const Vec3<S>& v)
    {
        const S x = v(0);
        const S y = v(1);
        const S z = v(2);
        Eigen::Matrix3<S> skew;
        skew <<
            S(0), -z, y,
            z, S(0), -x,
            -y, x, S(0);
        return skew;
    }

    template <typename S>
    Eigen::Matrix3<S> ToEigen(const SO3<S>& so3)
    {
        Eigen::Matrix3<S> R;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                R(i, j) = so3(i, j);
            }
        }
        return R;
    }

    struct OptimizeCameraIndex
    {
        int index2D{ -1 };
        int index3D{ -1 };
        int indexSE3{ -1 };
        int indexCam{ -1 };
    };

    struct ConstantParameter
    {
        enum Type
        {
            NONE,
            CAMERA_MODEL,
            CAMERA_POSE,
            DEVICE_POSE
        };

        Type type{ Type::NONE };
        int index{ 0 };
    };

    template <typename S>
    S OptimizeCameraIteration(
        const std::vector<OptimizeCameraIndex>& indices,
        const std::vector<Vec2<S>*>& points2D,
        const std::vector<Vec3<S>*>& points3D,
        std::vector<SE3<S>*>& devicesFromWorld,
        std::vector<SE3<S>*>& camsFromDevice,
        std::vector<ICameraModel<S>*>& cameras,
        S lambda = S(0),
        std::vector<ConstantParameter> constants = {})
    {
        using Eigen::Dynamic;
        using Eigen::RowMajor;
        using Eigen::Matrix;
        using Eigen::MatrixX;
        using Eigen::Vector;
        using Eigen::VectorX;
        using Eigen::Map;

        // If the constants are empty, we need at least
        // the first camera extrinsics to be constant.
        if (constants.empty())
        {
            constants.push_back({ ConstantParameter::CAMERA_POSE, 0 });
        }

        std::vector<bool> keepFrameConstant(devicesFromWorld.size(), false);
        std::vector<bool> keepCameraModelConstant(cameras.size(), false);
        std::vector<bool> keepCameraPoseConstant(camsFromDevice.size(), false);
        for (const auto& constant : constants)
        {
            if (constant.type == ConstantParameter::CAMERA_POSE)
            {
                keepCameraPoseConstant[constant.index] = true;
            }
            else if (constant.type == ConstantParameter::CAMERA_MODEL)
            {
                keepCameraModelConstant[constant.index] = true;
            }
            else if (constant.type == ConstantParameter::DEVICE_POSE)
            {
                keepFrameConstant[constant.index] = true;
            }
        }

        S totalSqResidual = S(0);

        int numCameraParams = 0;
        std::vector<int> cameraIndices;
        for (int i = 0; i < cameras.size(); ++i)
        {
            cameraIndices.push_back(numCameraParams);
            numCameraParams += cameras[i]->NumParams() + 6;
        }

        int numFrameParams = 0;
        std::vector<int> frameIndices;
        for (int i = 0; i < devicesFromWorld.size(); ++i)
        {
            frameIndices.push_back(numCameraParams + numFrameParams);
            numFrameParams += 6;
        }

        const int numTotalParams = numCameraParams + numFrameParams;
        Matrix<S, Dynamic, Dynamic, RowMajor> JTJ;
        JTJ.resize(numTotalParams, numTotalParams);
        JTJ.setZero();

        VectorX<S> g;
        g.resize(numTotalParams);
        g.setZero();

        // Convenience lambdas
        auto GetHessianBlock_CC = [&](int i) {
            const int n = cameras[i]->NumParams() + 6;
            return JTJ.block(cameraIndices[i], cameraIndices[i], n, n);
        };

        auto GetHessianBlock_CF = [&](int i, int j) {
            const int n = cameras[i]->NumParams() + 6;
            return JTJ.block(cameraIndices[i], frameIndices[j], n, 6);
        };

        auto GetHessianBlock_FC = [&](int i, int j) {
            const int n = cameras[j]->NumParams() + 6;
            return JTJ.block(frameIndices[i], cameraIndices[j], 6, n);
        };

        auto GetHessianBlock_FF = [&](int i) {
            return JTJ.block<6, 6>(frameIndices[i], frameIndices[i]);
        };

        auto GetGradientBlock_C = [&](int i) {
            return g.block(cameraIndices[i], 0, cameras[i]->NumParams() + 6, 1);
        };

        auto GetGradientBlock_F = [&](int i) {
            return g.block<6, 1>(frameIndices[i], 0);
        };

        // Single storage for common jacobians
        EigenList<Matrix<S, Dynamic, Dynamic, RowMajor>> J_uv_cams(cameras.size());
        EigenList<Matrix<S, Dynamic, Dynamic, RowMajor>> J_uv_models(cameras.size());
        for (int i = 0; i < cameras.size(); ++i)
        {
            J_uv_cams[i].resize(2, cameras[i]->NumParams() + 6);
            J_uv_models[i].resize(2, cameras[i]->NumParams());
        }

        Vec2<S> uv;
        Matrix<S, 2, 6> J_uv_frame;
        Matrix<S, 2, 3, RowMajor> J_uv_xyz;
        Matrix<S, 3, 6> J_xyz_frame;
        Matrix<S, 3, 6> J_xyz_extrinsics;

        auto FiniteDiff_Model = [&](auto& index)
        {
            const auto& point2D = *points2D[index.index2D];
            const auto& point3D = *points3D[index.index3D];
            const auto& deviceFromWorld = *devicesFromWorld[index.indexSE3];
            const auto& camFromDevice = *camsFromDevice[index.indexCam];
            const auto& camera = cameras[index.indexCam];

            MatrixX<S> J;
            J.resize(2, camera->NumParams());
            J.setZero();
            const auto eps = sqrt(std::numeric_limits<S>::epsilon());
            for (int i = 0; i < camera->NumParams(); ++i)
            {
                const auto point = camFromDevice * deviceFromWorld * point3D;
                auto cam_hi = camera->Clone();
                auto cam_lo = camera->Clone();

                cam_hi->Param(i) += eps;
                cam_lo->Param(i) -= eps;

                Vec2<S> uv_hi, uv_lo;
                cam_hi->Project(point, uv_hi);
                cam_lo->Project(point, uv_lo);

                Vec2<S> res = uv_hi - uv_lo;
                J(0, i) = res(0) / (2 * eps);
                J(1, i) = res(1) / (2 * eps);
            }

            return J;
        };
        
        auto FiniteDiff_Pose = [&](auto& index)
        {
            const auto& point2D = *points2D[index.index2D];
            const auto& point3D = *points3D[index.index3D];
            const auto& deviceFromWorld = *devicesFromWorld[index.indexSE3];
            const auto& camFromDevice = *camsFromDevice[index.indexCam];
            const auto& camera = cameras[index.indexCam];

            MatrixX<S> J;
            J.resize(2, 6);
            J.setZero();
            const auto eps = sqrt(std::numeric_limits<S>::epsilon());
            for (int i = 0; i < 6; ++i)
            {
                Vec2<S> uv_hi, uv_lo;

                Vec6<S> update_hi(S(0)), update_lo(S(0));
                update_hi(i) = eps;
                update_lo(i) = -eps;

                const auto point_hi = camFromDevice * SE3<S>::Exp(update_hi) * deviceFromWorld * point3D;
                const auto point_lo = camFromDevice * SE3<S>::Exp(update_lo) * deviceFromWorld * point3D;

                camera->Project(point_hi, uv_hi);
                camera->Project(point_lo, uv_lo);

                Vec2<S> res = uv_hi - uv_lo;
                J(0, i) = res(0) / (2 * eps);
                J(1, i) = res(1) / (2 * eps);
            }

            return J;
        };

        auto FiniteDiff_Extrinsics = [&](auto& index)
        {
            const auto& point2D = *points2D[index.index2D];
            const auto& point3D = *points3D[index.index3D];
            const auto& deviceFromWorld = *devicesFromWorld[index.indexSE3];
            const auto& camFromDevice = *camsFromDevice[index.indexCam];
            const auto& camera = cameras[index.indexCam];

            MatrixX<S> J;
            J.resize(2, 6);
            J.setZero();
            const auto eps = sqrt(std::numeric_limits<S>::epsilon());
            for (int i = 0; i < 6; ++i)
            {
                Vec2<S> uv_hi, uv_lo;

                Vec6<S> update_hi(S(0)), update_lo(S(0));
                update_hi(i) = eps;
                update_lo(i) = -eps;

                const auto point_hi = SE3<S>::Exp(update_hi) * camFromDevice * deviceFromWorld * point3D;
                const auto point_lo = SE3<S>::Exp(update_lo) * camFromDevice * deviceFromWorld * point3D;

                camera->Project(point_hi, uv_hi);
                camera->Project(point_lo, uv_lo);

                Vec2<S> res = uv_hi - uv_lo;
                J(0, i) = res(0) / (2 * eps);
                J(1, i) = res(1) / (2 * eps);
            }

            return J;
        };

        for (const auto& index : indices)
        {
            const auto& point2D = *points2D[index.index2D];
            const auto& point3D = *points3D[index.index3D];
            const auto& deviceFromWorld = *devicesFromWorld[index.indexSE3];
            const auto& camFromDevice = *camsFromDevice[index.indexCam];
            const auto& camera = cameras[index.indexCam];

            auto& J_uv_cam = J_uv_cams[index.indexCam];
            auto& J_uv_model_temp = J_uv_models[index.indexCam];
            auto& J_uv_model = J_uv_cam.block(0, 0, 2, camera->NumParams());
            auto& J_uv_extrinsics = J_uv_cam.block<2, 6>(0, camera->NumParams());

            const auto pointInDevice = deviceFromWorld * point3D;
            const auto pointInCam = camFromDevice * pointInDevice;

            camera->Project(pointInCam, uv, J_uv_model_temp.data(), J_uv_xyz.data());

            const auto residual = uv - point2D;
            totalSqResidual += residual.SquaredNorm();

            const Eigen::Matrix3<S> R_c = ToEigen(camFromDevice.R);
            J_xyz_frame.block<3, 3>(0, 0) = R_c * Skew(-pointInDevice);
            J_xyz_frame.block<3, 3>(0, 3) = R_c;
            J_uv_frame = J_uv_xyz * J_xyz_frame;

            J_xyz_extrinsics.block<3, 3>(0, 0) = Skew(-pointInCam);
            J_xyz_extrinsics.block<3, 3>(0, 3).setIdentity();
            J_uv_model = J_uv_model_temp;
            J_uv_extrinsics = J_uv_xyz * J_xyz_extrinsics;

#if 0
            const auto J_uv_model_finite = FiniteDiff_Model(index);
            const auto delta_model = J_uv_model - J_uv_model_finite;
            bool print_line = false;
            const auto eps = std::numeric_limits<S>::epsilon();
            if (delta_model.squaredNorm() > eps)
            {
                std::cout << "J model delta" << std::endl;
                std::cout << J_uv_model << std::endl;
                std::cout << J_uv_model_finite << std::endl;
                std::cout << delta_model << std::endl;
                print_line = true;
            }

            const auto J_uv_extrinsics_finite = FiniteDiff_Extrinsics(index);
            const auto delta_extrinsics = J_uv_extrinsics - J_uv_extrinsics_finite;
            if (delta_extrinsics.squaredNorm() > eps)
            {
                std::cout << "J extrinsics delta" << std::endl;
                std::cout << J_uv_extrinsics << std::endl;
                std::cout << J_uv_extrinsics_finite << std::endl;
                std::cout << delta_extrinsics << std::endl;
                print_line = true;
            }

            const auto J_uv_frame_finite = FiniteDiff_Pose(index);
            const auto delta_frame = J_uv_frame - J_uv_frame_finite;
            if (delta_frame.squaredNorm() > eps)
            {
                std::cout << "J frame delta" << std::endl;
                std::cout << J_uv_frame << std::endl;
                std::cout << J_uv_frame_finite << std::endl;
                std::cout << delta_frame << std::endl;
                print_line = true;
            }
            if (print_line)
            {
                for (int i = 0; i < 3; ++i)
                {
                    std::cout << point3D(i) << ' ';
                }
                std::cout << std::endl;
                std::cout << "============================================" << std::endl;
            }
#endif

            const auto& camIdx = index.indexCam;
            const auto& se3Idx = index.indexSE3;
            GetHessianBlock_CC(camIdx) += J_uv_cam.transpose() * J_uv_cam;
            GetHessianBlock_FF(se3Idx) += J_uv_frame.transpose() * J_uv_frame;
            GetHessianBlock_CF(camIdx, se3Idx) += J_uv_cam.transpose() * J_uv_frame;
            GetHessianBlock_FC(se3Idx, camIdx) += J_uv_frame.transpose() * J_uv_cam;

            const Vector<S, 2> eResidual(residual(0), residual(1));
            GetGradientBlock_C(camIdx) += J_uv_cam.transpose() * eResidual;
            GetGradientBlock_F(se3Idx) += J_uv_frame.transpose() * eResidual;
        }

        // Dampen hessian
        for (int i = 0; i < JTJ.rows(); ++i)
        {
            JTJ(i, i) += lambda;
        }

        // Negate g now for simplicity
        g = -g;

        // [ A  B] [x0] = [y0]
        // [B^T C] [x1]   [y1]
        //
        //  A  * x0 + B * x1 = y0
        // B^T * x0 + C * x1 = y1
        //
        // x0 = (A - B * C^-1 * B^T)^-1 * (y0 - B * C^-1 * y1)
        // x1 = C^-1 * (y1 - B^T * x0)
        const MatrixX<S>& A = JTJ.block(0, 0, numCameraParams, numCameraParams);
        const MatrixX<S>& B = JTJ.block(0, numCameraParams, numCameraParams, numFrameParams);
        const MatrixX<S>& C = JTJ.block(numCameraParams, numCameraParams, numFrameParams, numFrameParams);

        const VectorX<S>& y0 = g.block(0, 0, numCameraParams, 1);
        const VectorX<S>& y1 = g.block(numCameraParams, 0, numFrameParams, 1);

        // Since C is block diagonal, its inverse is
        // also block diagonal with the inverse of each block.
        MatrixX<S> Cinv;
        Cinv.resize(numFrameParams, numFrameParams);
        Cinv.setZero();

        const Eigen::Matrix<S, 6, 6> I6 = Matrix<S, 6, 6>::Identity();
        for (int i = 0; i < devicesFromWorld.size(); ++i)
        {
            const int j = frameIndices[i] - numCameraParams;
            auto& C_block = C.block<6, 6>(j, j);
            auto& Cinv_block = Cinv.block<6, 6>(j, j);

            Cinv_block = C_block.ldlt().solve(I6);
        }

        // Compute x0 using full inverse of C
        const MatrixX<S> M = A - B * Cinv * B.transpose();
        const VectorX<S> b = y0 - B * Cinv * y1;
        const VectorX<S> x0 = M.ldlt().solve(b);

        // Solve for x1 by using LDLT solve per frame, given x0
        VectorX<S> x1;
        x1.resize(numFrameParams);
        x1.setZero();
        const VectorX<S> y1_prime = y1 - B.transpose() * x0;
        for (int i = 0; i < devicesFromWorld.size(); ++i)
        {
            const int j = frameIndices[i] - numCameraParams;
            const auto& y1p_i = y1_prime.block<6, 1>(j, 0);

            auto& x1_i = x1.block<6, 1>(j, 0);
            x1_i = GetHessianBlock_FF(i).ldlt().solve(y1p_i);
        }

        // Update all of the camera parameters
        for (int i = 0; i < cameras.size(); ++i)
        {
            auto& camera = cameras[i];
            auto& extrinsics = *camsFromDevice[i];
            int index = cameraIndices[i];
            if (!keepCameraModelConstant[i])
            {
                for (int j = 0; j < camera->NumParams(); ++j)
                {
                    camera->Param(j) += x0(index + j);
                }
            }

            if (!keepCameraPoseConstant[i])
            {
                index += camera->NumParams();
                const Eigen::Vector<S, 6> dt = x0.block<6, 1>(index, 0);
                const auto dT = SE3<S>::Exp(
                    { dt(0), dt(1), dt(2), dt(3), dt(4), dt(5) });
                extrinsics = dT * extrinsics;
            }
        }

        // Update all of the pose parameters
        for (int i = 0; i < devicesFromWorld.size(); ++i)
        {
            if (keepFrameConstant[i])
            {
                continue;
            }

            auto& pose = *devicesFromWorld[i];
            const int index = frameIndices[i] - numCameraParams;

            const Eigen::Vector<S, 6> dt = x1.block<6, 1>(index, 0);
            const auto dT = SE3<S>::Exp(
                { dt(0), dt(1), dt(2), dt(3), dt(4), dt(5) });

            pose = dT * pose;
        }

        return totalSqResidual;
    }
}
