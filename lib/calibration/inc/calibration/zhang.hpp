#pragma once

#include <vector>
#include <Eigen/SVD>
#include <lie/se3.hpp>
#include <calibration/pinholemodel.hpp>

namespace sight
{
    namespace
    {
        template <typename S>
        Eigen::Vector<S, 6> GetZhangConstraintVector(
            const Eigen::Vector3<S>& hi,
            const Eigen::Vector3<S>& hj)
        {
            Eigen::Vector<S, 6> v;
            v <<
                hi(0) * hj(0),
                hi(0) * hj(1) + hi(1) * hj(0),
                hi(1) * hj(1),
                hi(2) * hj(0) + hi(0) * hj(2),
                hi(2) * hj(1) + hi(1) * hj(2),
                hi(2) * hj(2);
            return v;
        }

        template <typename S>
        PinholeModel<S> IntrinsicsFromBVector(const Eigen::Vector<S, 6>& B)
        {
            PinholeModel<S> K;

            const S B11 = B(0);
            const S B12 = B(1);
            const S B13 = B(2);
            const S B22 = B(3);
            const S B23 = B(4);
            const S B33 = B(5);

            const S cy = (B12 *B13 - B11 * B23) / (B11 * B22 - B12 * B12);
            const S lambda = B33 - (B13 * B13 + cy * (B12 * B13 - B11 * B23)) / B11;
            const S fx = sqrt(lambda / B11);
            const S fy = sqrt((lambda * B11) / (B11 * B22 - B12 * B12));
            const S skew = -B12 * fx * fx * fy / lambda;
            const S cx = skew * cy / fy - B13 * fx * fx / lambda;

            K.param(PinholeModel<S>::FX) = fx;
            K.param(PinholeModel<S>::FY) = fy;
            K.param(PinholeModel<S>::CX) = cx;
            K.param(PinholeModel<S>::CY) = cy;
            K.param(PinholeModel<S>::SKEW) = skew;

            return K;
        }
    }

    template <typename S>
    Eigen::Matrix3<S> ClosestFrobeniusRotation(const Eigen::Matrix3<S>& R_noisy)
    {
        const Eigen::JacobiSVD<Eigen::MatrixX<S>> svd(R_noisy, Eigen::ComputeThinU | Eigen::ComputeThinV);
        return svd.computeU() * svd.computeV().transpose();
    }

    template <typename S>
    void PoseFromMetricHomography(
        const Eigen::Matrix3<S>& H,
        SE3<S>& dstFromPlaneZ0)
    {
        const Eigen::Vector3<S> sR1 = H.col(0);
        const Eigen::Vector3<S> sR2 = H.col(1);
        const Eigen::Vector3<S> st  = H.col(2);

        const S s = (sR1.norm() + sR2.norm()) / S(2);
        const S sinv = S(1) / s;

        Eigen::Matrix3<S> R_noisy;
        R_noisy.col(0) = sR1.normalized();
        R_noisy.col(1) = sR2.normalized();
        R_noisy.col(2) = (R_noisy.col(0).cross(R_noisy.col(1))).normalized();

        const Eigen::Matrix3<S> R = ClosestFrobeniusRotation(R_noisy);
        const Eigen::Vector3<S> t = st * sinv;

        // Copy into the SE(3) group element
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                dstFromPlaneZ0.R(i, j) = R(i, j);
            }
            dstFromPlaneZ0.t(i) = t(i);
        }
    }

    template <typename S>
    void PoseFromImageHomography(
        const Eigen::Matrix3<S>& H,
        const PinholeModel<S>& intrinsics,
        SE3<S>& dstFromPlaneZ0)
    {
        const S fx = intrinsics.param(PinholeModel<S>::FX);
        const S fy = intrinsics.param(PinholeModel<S>::FY);
        const S cx = intrinsics.param(PinholeModel<S>::CX);
        const S cy = intrinsics.param(PinholeModel<S>::CY);
        const S skew = intrinsics.param(PinholeModel<S>::SKEW);

        Eigen::Matrix3<S> K;
        K <<
            fx, skew, cx,
            S(0), fy, cy,
            S(0), s(0), S(1);
        
        const Eigen::Matrix3<S> Kinv = K.inverse();
        const Eigen::Matrix3<S> sR1R2t = Kinv * H;

        return PoseFromMetricHomography(sR1R2t);
    }

    template <typename S>
    PinholeModel<S> ZhangIntrinsicsFromHomographies(
        const std::vector<Eigen::Matrix3<S>>& Hs)
    {
        // To use Zhang's method of extracting pose, [R t], from
        // a homography, H, we assume that all of the points that
        // are mapped by the homography come from a planar surface.
        //
        // Thus, in the tracked planar surface's coordinate frame,
        // its points lie on the z=0 plane -- or all points are
        // of the form [x, y, 0].
        //
        // Thus, with a linear calibration, we have:
        //
        //   uv{cam} = 
        //   uv{cam} = s * K * [ r1 | r2 | t ] * [x, y, 0]^T
        //
        // Where r1 and r2 are the first columns of the rotation
        // matrix and t is the translation. Note, that homographies
        // are defined only up to a scale. So there's a scalar, s,
        // leftover.
        //
        // We have constraints on r1 and r2 -- they must be orthogonal.
        // And they also must have the same norm (since they're both
        // unit 1 and scaling it doesn't change that).
        //
        // So,
        //    (K^-1 * H.col(0))^T * (K^-1 * H.col(1)) = 0
        //    (K^-1 * H.col(0))^T(K^-1 * H.col(0)) = (K^-1 * H.col(1))^T(K^-1 * H.col(1))
        //
        // All of the operations share the same matrix:
        //
        //    A = K^-T * K^-1
        //
        // Thus, rewriting our constraints:
        //
        //    H.col(0)   * A * H.col(1) = 0
        //    H.col(0)^T * A * H.col(0) = H.col(1)^T * A * H.col(1)
        //
        // We get 2 constraints for A, which is a symmetric 3x3 matrix.
        //
        // There are only 6 unique values in a 3x3 symmetric matrix, so
        // we require 3 homographies to estimate A.
        const int n = Hs.size();
        Eigen::MatrixX<S> M;
        M.resize(n, 6);

        for (int i = 0; i < n; ++i)
        {
            const auto& H = Hs[i];
            M.row(i * 2 + 0) = GetZhangConstraintVector(H.col(0), H.col(1));
            M.row(i * 2 + 1) =
                GetZhangConstraintVector(H.col(0), H.col(0)) -
                GetZhangConstraintVector(H.col(1), H.col(1));
        }

        const Eigen::JacobiSVD<Eigen::MatrixX<S>> svd(M, Eigen::ComputeFullV);
        const Eigen::Vector<S, 6> B = svd.computeV().col(5);

        return IntrinsicsFromBVector(B);
    }
}