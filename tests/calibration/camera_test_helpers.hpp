#pragma once

#include <memory>
#include <iostream>

#include <catch2/catch.hpp>
#include <calibration/cameramodel.hpp>
#include <linear/vec.hpp>

template <typename S>
void TestCamera(const sight::CameraModel<S>& cam)
{
    using namespace sight;
    
    // Create a cube-grid of points in the range [-.5, .5] x [-.5, .5] x [.25, 1.25]
    const int N = 20;
    const S step = S(1) / (N - 1);
    
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int k = 0; k < N; ++k)
            {
                const S x = S(-.5) + i * step;
                const S y = S(-.5) + j * step;
                const S z = S(.25) + k * step;

                std::vector<S> Juv_xyz(2 * 3);
                std::vector<S> Juv_xzyz(2 * 2);
                std::vector<S> Juv_params(2 * cam.NumParams());

                S u, v;
                cam.Project(x, y, z, u, v, Juv_params.data(), Juv_xyz.data());
                cam.JxyzToJxzyz(z, Juv_xyz.data(), Juv_xzyz.data());

                const S eps = sqrt(std::numeric_limits<S>::epsilon());

                // Test d{u,v} / d{params}
                {
                    std::vector<S> Juv_params_finite(Juv_params.size());

                    // Perturb each parameter
                    std::unique_ptr<CameraModel<S>> cam_lo, cam_hi;
                    for (int idx = 0; idx < cam.NumParams(); ++idx)
                    {
                        // Re-copy the current camera parameters over
                        cam_lo = cam.Clone();
                        cam_hi = cam.Clone();

                        cam_lo->Param(idx) -= eps;
                        cam_hi->Param(idx) += eps;

                        Vec2<S> uv_lo, uv_hi;
                        cam_lo->Project(x, y, z, uv_lo[0], uv_lo[1]);
                        cam_hi->Project(x, y, z, uv_hi[0], uv_hi[1]);

                        Juv_params_finite[idx + 0] = (uv_hi[0] - uv_lo[0]) / (2 * eps);
                        Juv_params_finite[idx + cam.NumParams()] = (uv_hi[1] - uv_lo[1]) / (2 * eps);
                    }

                    for (int idx = 0; idx < Juv_params.size(); ++idx)
                    {
                        REQUIRE(Juv_params[idx] == Approx(Juv_params_finite[idx]).margin(1e-4));
                    }
                }

                // Test d{u,v}/d{x, y, z}
                {
                    std::vector<S> Juv_xyz_finite(2 * 3);
                    Vec3<S> xyz_lo, xyz_hi;
                    for (int idx = 0; idx < 3; ++idx)
                    {
                        xyz_lo = {x, y, z};
                        xyz_hi = {x, y, z};

                        xyz_lo[idx] -= eps;
                        xyz_hi[idx] += eps;

                        Vec2<S> uv_lo, uv_hi;
                        cam.Project(xyz_lo[0], xyz_lo[1], xyz_lo[2], uv_lo[0], uv_lo[1]);
                        cam.Project(xyz_hi[0], xyz_hi[1], xyz_hi[2], uv_hi[0], uv_hi[1]);

                        Juv_xyz_finite[idx + 0] = (uv_hi[0] - uv_lo[0]) / (2 * eps);
                        Juv_xyz_finite[idx + 3] = (uv_hi[1] - uv_lo[1]) / (2 * eps);
                    }

                    for (int idx = 0; idx < 6; ++idx)
                    {
                        REQUIRE(Juv_xyz[idx] == Approx(Juv_xyz_finite[idx]).margin(1e-4));
                    }
                }

                // Test d{u,v}/d{x/z,y/z}
                {
                    std::vector<S> Juv_xzyz_finite(2 * 2);
                    Vec2<S> xzyz_lo, xzyz_hi;
                    for (int idx = 0; idx < 2; ++idx)
                    {
                        xzyz_lo = {x / z, y / z};
                        xzyz_hi = {x / z, y / z};

                        xzyz_lo[idx] -= eps;
                        xzyz_hi[idx] += eps;

                        Vec2<S> uv_lo, uv_hi;
                        cam.Project(xzyz_lo[0], xzyz_lo[1], S(1), uv_lo[0], uv_lo[1]);
                        cam.Project(xzyz_hi[0], xzyz_hi[1], S(1), uv_hi[0], uv_hi[1]);

                        Juv_xzyz_finite[idx + 0] = (uv_hi[0] - uv_lo[0]) / (2 * eps);
                        Juv_xzyz_finite[idx + 2] = (uv_hi[1] - uv_lo[1]) / (2 * eps);
                    }

                    for (int idx = 0; idx < 4; ++idx)
                    {
                        REQUIRE(Juv_xzyz[idx] == Approx(Juv_xzyz_finite[idx]).margin(1e-4));
                    }
                }
            
                // Test unprojection
                {
                    S xz, yz;
                    const bool success = cam.Unproject(u, v, xz, yz);

                    S up, vp;
                    cam.Project(xz, yz, S(1), up, vp);
                    
                    REQUIRE(success);
                    REQUIRE(xz == Approx(x / z).margin(eps));
                    REQUIRE(yz == Approx(y / z).margin(eps));

                    REQUIRE(up == Approx(u).margin(eps));
                    REQUIRE(vp == Approx(v).margin(eps));
                }
            }
        }
    }
}