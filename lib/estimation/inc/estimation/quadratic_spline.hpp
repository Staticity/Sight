#pragma once

#include <numeric>
#include <algorithm>

#include <estimation/quadratic.hpp>

namespace sight
{

    template <typename S>
    class QuadraticSpline
    {
    public:

        QuadraticSpline()
            : xs()
            , ys()
            , indices()
            , quads()
            , built(false)
        {}

        void AddPoint(S x, S y)
        {
            xs.push_back(x);
            ys.push_back(y);
            built = false;
        }

        bool IsBuilt() const
        {
            return built;
        }
        
        void BuildSpline()
        {
            // Reset previously built quadratics
            built = false;
            quads.clear();
            indices.clear();

            // Insufficient points?
            if (xs.size() < 3)
            {
                return;
            }

            indices.resize(xs.size());
            quads.resize(xs.size() - 2);

            // Since points could have been added in arbitrary order, 
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](int i, int j) { return xs[i] < xs[j]; });

            for (int idx = 1; idx < quads.size() - 1; ++idx)
            {
                const int i = indices[idx - 1];
                const int j = indices[idx];
                const int k = indices[idx + 1];
                quads[idx] = Quadratic<S>(xs[i], ys[i], xs[j], ys[j], xs[k], ys[k]);
            }

            built = true;
        }

        void GetExtrema(
            std::vector<S>& outXs,
            std::vector<S>& outYs,
            std::vector<bool>& isMax)
        {
            if (!built)
            {
                BuildSpline();
            }

            for (int idx = 0; idx < quads.size(); ++idx)
            {
                const auto& quad = quads[idx];
                const S xMin = xs[indices[idx]];
                const S xMax = xs[indices[idx + 2]];

                const S yMin = ys[indices[idx]];
                const S yMax = ys[indices[idx + 2]];

                S x;
                if (quad.GetExtremum(x))
                {
                    // Check if the extremum
                    if (x >= xMin && x <= xMax && (outXs.empty() || outXs.back() != x))
                    {
                        outXs.push_back(x);
                        outYs.push_back(quad.Eval(x));
                        isMax.push_back(quad.HasMaximum());
                    }
                }
            }
        }

        bool built;        

        std::vector<S> xs;
        std::vector<S> ys;

        std::vector<int> indices;
        std::vector<Quadratic<S>> quads;
    };
}