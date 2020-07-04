#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <image/image_gen.hpp>
#include <linear/vec.hpp>

namespace sight
{

    template <typename S>
    struct Flow
    {
        bool valid;
        int radius;
        Vec2<S> pos;
        Vec2<S> velocity;
        Eigen::Matrix3<S> warp;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
    
    template <typename T, typename S>
    void OpticalFlow_KLT(
        const Image<T>& image0,
        const Image<T>& image1,
        std::vector<Flow<S>>& flows,
        const int pxOffset = 5,
        const int localRadius = 5,
        const int maxDisplacement = 21,
        const int maxRefinementIterations = 10,
        const S minEigenValue = S(.01))
    {
        // Only supporting grayscale images for now.
        assert(image0.c == 1);
        assert(image1.c == 1);
        if (image0.c != 1 || image1.c != 1)
        {
            return;
        }

        const Image<T> view0 = PadView(image0, localRadius + 1);

        // Not necessary to compute the entirety of the image's gradient.
        // Maybe it's faster to compute each region itself (although, that
        // would be slower for dense flows).
        const int padding = std::max({localRadius, maxDisplacement}) + 1;
        const Image<int> dx = PadView(Sobel3x3_Horizontal<int, T>(view0), padding);
        const Image<int> dy = PadView(Sobel3x3_Vertical<int, T>(view0), padding);
        const Image<T> view1 = PadView(image1, padding);

        std::vector<S> gKernel = GaussianKernel<S>(localRadius, 2.0);
        std::vector<S> weights(gKernel.size() * gKernel.size());
        for (int i = 0, k = 0; i < gKernel.size(); ++i)
        {
            for (int j = 0; j < gKernel.size(); ++j, ++k)
            {
                weights[k] = S(1); //gKernel[i] * gKernel[j];
                weights[k] *= weights[k];
            }
        }

        const int numExpectedFlows =
            ((image0.w + pxOffset - 1) / pxOffset) * ((image0.h + pxOffset - 1) / pxOffset);
        flows.reserve(numExpectedFlows);

        const int dispThreshSq = maxDisplacement * maxDisplacement;
        const int sx = pxOffset / 2;
        const int sy = sx;
        for (int y = sy; y < view1.h; y += pxOffset)
        {
            for (int x = sx; x < view1.w; x += pxOffset)
            {
                // We'll use for efficiency's sake later (maybe)
                // const T* src = view0.at(y, x);

                // We'll solve for a translational transformation initially,
                // just so we can test a simpler implementation.
                
                // W is parameterized by w, where
                //
                //     w = {tx, ty}
                //
                // W will 'warp' a point p = {x, y}
                //
                //     W(p;w) = [x + tx, y + ty]
                //
                // And its derivative w.r.t p is:
                //
                //     W'(p;w) = [1 0]
                //               [0 1]
                //

                Flow<S> flow;
                flow.valid = true;
                flow.radius = localRadius;
                flow.pos = { S(x), S(y) };
                flow.velocity.Fill(S(0));
                flow.warp = Eigen::Matrix3<S>::Identity();

                // Initialize the translation to 0
                Vec2<S> w = { S(0), S(0) };

                for (int i = 0; i < maxRefinementIterations; ++i)
                {
                    const S centerDispSq = w.SquaredNorm();
                    if (centerDispSq > dispThreshSq)
                    {
                        flow.valid = false;
                        break;
                    }

                    // Should be NxN where N is dimension of the warping function
                    Eigen::Matrix<S, 2, 2> JTJ;
                    JTJ.setZero();

                    // Should be 1xN where N is dimension of the warping function
                    Eigen::Vector<S, 2> g;
                    g.setZero();

                    S totalResSq = S(0);

                    // For all pixels in within [-r, r] x [-r, r] centered
                    // at this pixel..
                    for (int oy = -localRadius, k = 0; oy <= localRadius; ++oy)
                    {
                        for (int ox = -localRadius; ox <= localRadius; ++ox, ++k)
                        {
                            // Compute the warped point
                            const Vec2<int> p = { x + ox, y + oy };
                            const Vec2<S> Wp = { p[0] + w[0], p[1] + w[1] };

                            // Compute the warped point's intensity
                            const S pred = BilinearInterpolate<S, T, S>(view1, Wp[0], Wp[1]);
                            const S obs = view0(p[1], p[0]);

                            // Compute Jacobian of this residual
                            // Since, we're using translation, it's
                            // just the Jacobian of the original image.
                            Eigen::Matrix<S, 1, 2> J;
                            J << S(dx(p[1], p[0])), S(dy(p[1], p[0]));

                            const S weight = weights[k];
                            const S residual = weight * (pred - obs);
                            totalResSq += residual * residual;

                            JTJ += J.transpose() * J;
                            g += J.transpose() * residual;
                        }
                    }

                    Eigen::Vector2<std::complex<S>> eigenvalues = JTJ.eigenvalues();
                    if (std::min(eigenvalues[0].real(), eigenvalues[1].real()) <= minEigenValue)
                    {
                        flow.valid = false;
                        break;
                    }

                    // std::cout << i << " " << totalResSq / pow(localRadius * 2 + 1, 2) << std::endl;

                    const Eigen::Matrix<S, 2, 1> wDelta = JTJ.colPivHouseholderQr().solve(-g);
                    w[0] += wDelta(0);
                    w[1] += wDelta(1);

                    // convergence
                    if (wDelta.squaredNorm() <= std::numeric_limits<S>::epsilon())
                    {
                        break;
                    }
                }

                flow.velocity = { w[0], w[1] };
                // if (flow.valid)
                // {
                //     // flow.warp = ...
                // }

                flows.push_back(flow);
            }
        }
    }

}