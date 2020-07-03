#include <iostream>
#include <map>
#include <fstream>
#include <filesystem>

#include <Eigen/Eigen>
#include <lie/so3.hpp>
#include <lie/lie_ops.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace sight;

// template <typename S>
// S FindSmallest()
// {
//     S v = S(1);
//     while (v) v /= S(2);
//     return v;
// }

// N evenly spaced values from (lo, hi]
template <typename S>
std::vector<S> Linspace(S lo, S hi, int N)
{
    const S inc = (hi - lo) / N;
    std::vector<S> v(N);
    for (int i = 0; i < N; ++i)
    {
        v[i] = lo;
        lo += inc;
    }
    return v;
}

template <typename S>
S FindOptimalThreshold(
    const S minThresh,
    const S maxThresh,
    const int nThresholds,
    const S minTheta,
    const S maxTheta,
    const int nThetas,
    bool debug = true)
{
    const std::string precision = std::string(typeid(S).name());

    namespace fs = std::filesystem;
    const fs::path curr = __FILE__;
    const fs::path out_dir = curr.parent_path() / "outputs" / precision;

    const fs::path output_path = out_dir / "result.csv";
    std::ofstream out(output_path.string());
    out.precision(std::numeric_limits<double>::max_digits10);
    out << "Threshold,TotalErr\n";

    const fs::path debug_path = out_dir / "verbose.csv";
    std::ofstream debug_out(debug_path.string());

    // We'll output each individual residual for debugging
    debug_out.precision(std::numeric_limits<double>::max_digits10);
    debug_out << "threshIdx,theta,id,err\n";

    const std::vector<S> threshes = Linspace<S>(minThresh, maxThresh, nThresholds);
    const std::vector<S> thetas = Linspace<S>(minTheta, maxTheta, nThetas);

    S bestThresh = threshes[0];
    double bestError = std::numeric_limits<double>::infinity();

    // Hard-coded for now:
    //
    // Store a set of angle axes to exponentiate with
    // the brute-force method and SO3::Exp method.
    const int nVecsPerTheta = 1;//3;
    std::map<int, std::vector<Vec3<S>>> angleAxes;
    for (int i = 0; i < thetas.size(); ++i)
    {
        const S theta = thetas[i];
        angleAxes[i].push_back({theta, 0, 0});
        angleAxes[i].push_back({0, theta, 0});
        angleAxes[i].push_back({0, 0, theta});
    }

    // Compute the ground truth results for each angle axis.
    std::map<int, std::vector<Eigen::Matrix3<double>>> expected;
    for (const auto& x : angleAxes)
    {
        const int idx = x.first;
        for (const auto& aa : x.second)
        {
            Eigen::Matrix3<double> A;
            A <<
                     0, -aa(2),  aa(1),
                 aa(2),      0, -aa(0),
                -aa(1),  aa(0),      0;
            
            expected[idx].push_back(Exp<double, 3>(A));
        }
    }

    // For each threshold...
    for (int threshIdx = 0; threshIdx < threshes.size(); ++threshIdx)
    {
        const S thresh = threshes[threshIdx];
        double totalErr = S(0);
        int numResiduals = 0;

        // Consider every possible theta
        for (int thetaIdx = 0; thetaIdx < thetas.size(); ++thetaIdx)
        {
            const S theta = thetas[thetaIdx];
            const auto& axes = angleAxes[thetaIdx];
            const auto& actuals = expected[thetaIdx];

            // For each angle axes with magnitude theta
            for (int i = 0; i < axes.size(); ++i)
            {
                const Vec3<S>& aa = axes[i];
                SO3<S> R = SO3<S>::Exp(aa, thresh);

                Eigen::Matrix3<double> estimated;
                estimated <<
                    R.R[0], R.R[1], R.R[2],
                    R.R[3], R.R[4], R.R[5],
                    R.R[6], R.R[7], R.R[8];
                const Eigen::Matrix3<double>& actual = actuals[i];

                const double trace = (estimated.transpose() * actual).trace();
                const double angularErr = acos(std::clamp((trace - 1) / 2.0, -1.0, 1.0));

                totalErr += angularErr;
                ++numResiduals;

                if (debug)
                {
                    debug_out << threshIdx << ',' << theta << ',' << i << ',' << angularErr << '\n';
                }
            }
        }
    
        out << thresh << "," << totalErr / numResiduals << '\n';

        if (totalErr < bestError)
        {
            bestError = totalErr;
            bestThresh = thresh;
        }
    }

    return bestThresh;
}

int main()
{   
    const bool debug = true;
    
    const float bestF =
        FindOptimalThreshold<float>(
            std::numeric_limits<float>::epsilon(),
            .002f,
            1000,
            0.0f,
            .001f,
            1000,
            debug);
    std::cout << "Best Floating Point Threshold: " << bestF << std::endl;

    const double bestD =
        FindOptimalThreshold<double>(
            std::numeric_limits<double>::epsilon(),
            1e-7,
            1000,
            0.0,
            1e-8,
            1000,
            debug);

    std::cout << "Best Double Floating Point Threshold: " << bestD << std::endl;
}