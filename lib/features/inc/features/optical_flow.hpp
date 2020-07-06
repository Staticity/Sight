#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <image/image_gen.hpp>
#include <linear/vec.hpp>
#include <utils/timer.hpp>

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
    
    // TODO: Add affine warping?
    template <typename T, typename S>
    void OpticalFlow_KLT(
        const Image<T>& image0,
        const Image<T>& image1,
        std::vector<Flow<S>>& flows,
        const bool blurImage = true,
        const int pxOffset = 5,
        const int localRadius = 5,
        const int maxDisplacement = 21,
        const int maxRefinementIterations = 10,
        const S minEigenValue = S(.01),
        const S updateEpsilon = std::numeric_limits<S>::epsilon())
    {
        // Only supporting grayscale images for now.
        assert(image0.c == 1);
        assert(image1.c == 1);
        if (image0.c != 1 || image1.c != 1)
        {
            return;
        }

        // Radius for gaussian blurring
        const int gaussR = 5;

        Image<T> view0;
        if (blurImage)
        {
            view0 = PadView(GaussianBlur<T, T>(PadView(image0, gaussR), 1.0f, gaussR), localRadius + 1);
        }
        else
        {
            view0 = PadView(image0, localRadius + 1);
        }

        // Not necessary to compute the entirety of the image's gradient.
        // Maybe it's faster to compute each region itself (although, that
        // would be slower for dense flows).
        using TT = std::common_type<int, T>::type;
        const int padding = localRadius + maxDisplacement + 1;
        const Image<TT> Idx = PadView(Sobel3x3_Horizontal<TT, T>(view0), padding);
        const Image<TT> Idy = PadView(Sobel3x3_Vertical<TT, T>(view0), padding);
        Image<T> view1;
        if (blurImage)
        {
            view1 = PadView(GaussianBlur<T, T>(PadView(image1, gaussR), 1.0f, gaussR), padding);
        }
        else
        {
            view1 = PadView(image1, padding);            
        }

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

        // To be re-used often..
        //
        // Should be NxN where N is dimension of the warping function
        // Should be 1xN where N is dimension of the warping function
        Vec<S, 2 * 2> JTJ;
        Vec2<S> g;

        const int dispThreshSq = maxDisplacement * maxDisplacement;
        const int sx = pxOffset / 2;
        const int sy = sx;
        for (int y = sy; y < view1.h; y += pxOffset)
        {
            for (int x = sx; x < view1.w; x += pxOffset)
            {
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
                flow.pos[0] = S(x);
                flow.pos[1] = S(y);
                flow.velocity[0] = S(0);
                flow.velocity[1] = S(0);
                // flow.warp = Eigen::Matrix3<S>::Identity();

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

                    // Reset our Hessian and gradient for each iteration
                    JTJ.Zero();
                    g.Zero();

                    // We'll use a pointer to the row for efficiency's sake
                    const T* srcRow = view0.at(y - localRadius, x);

                    // For all pixels in within [-r, r] x [-r, r] centered
                    // at this pixel..
                    for (int oy = -localRadius, k = 0; oy <= localRadius; ++oy, srcRow += view0.row_step)
                    {
                        for (int ox = -localRadius; ox <= localRadius; ++ox, ++k)
                        {
                            // Compute the warped point
                            // const Vec2<int> p = { x + ox, y + oy };
                            // const Vec2<S> Wp = { p[0] + w[0], p[1] + w[1] };
                            const int px = x + ox;
                            const int py = y + oy;
                            const S Wpx = px + w[0];
                            const S Wpy = py + w[1];

                            // Compute the warped point's intensity
                            const S pred = BilinearInterpolate<S, T, S>(view1, Wpx, Wpy);
                            const S obs = srcRow[ox];

                            // Compute Jacobian of this residual
                            // Since, we're using translation, it's
                            // just the Jacobian of the original image.
                            const TT dx = Idx(py, px);
                            const TT dy = Idy(py, px);
                            const TT dxdy = dx * dy;

                            const S weight = weights[k];
                            const S residual = weight * (pred - obs);

                            // JTJ += J.Transpose() * J;
                            JTJ[0] += dx * dx;
                            JTJ[1] += dxdy;
                            JTJ[2] += dxdy;
                            JTJ[3] += dy * dy;

                            // g += J.transpose() * residual;
                            g[0] += dx * residual;
                            g[1] += dy * residual;
                        }
                    }

                    Eigen::Matrix2<S> JTJe;
                    JTJe << JTJ[0], JTJ[1], JTJ[2], JTJ[3];
                    const Eigen::Vector2<S> ge(g[0], g[1]);

                    Eigen::Vector2<std::complex<S>> eigenvalues = JTJe.eigenvalues();
                    if (std::min(eigenvalues[0].real(), eigenvalues[1].real()) <= minEigenValue)
                    {
                        flow.valid = false;
                        break;
                    }

                    const Eigen::Matrix<S, 2, 1> wDelta = JTJe.colPivHouseholderQr().solve(-ge);
                    w[0] += wDelta(0);
                    w[1] += wDelta(1);

                    // convergence
                    if (wDelta.squaredNorm() <= updateEpsilon)
                    {
                        break;
                    }
                }

                flow.velocity[0] = w[0];
                flow.velocity[1] = w[1];
                flows.push_back(flow);
            }
        }
    }
}