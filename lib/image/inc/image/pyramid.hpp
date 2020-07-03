#pragma once

#include <image/image.hpp>
#include <image/image_ops.hpp>
#include <opencv2/imgproc.hpp>

namespace sight
{

    template <typename T>
    class Pyramid
    {
    public:

        Pyramid(const Image<T>& image, size_t maxLevels, const float scale)
        {
            // Currently only supporting grayscale pyramids.
            const size_t levels = std::min(maxLevels, ComputeMaxLevels(scale, image.w, image.h));

            images.resize(levels);
            scalesX.resize(levels, 1.f);
            scalesY.resize(levels, 1.f);

            images[0] = image;
            for (int i = 1; i < levels; ++i)
            {
                images[i] = Resize<T>(images[i - 1], scale, scale);

                // Due to truncation, the true scales may differ
                // slightly, so we recover that information.
                const float dsx = float(images[i - 1].w) / images[i].w;
                const float dsy = float(images[i - 1].h) / images[i].h;
                scalesX[i] = scalesX[i - 1] * dsx;
                scalesY[i] = scalesY[i - 1] * dsy;
            }
        }

        inline Image<T>& GetLevel(int i) { return images[i]; }
        inline const Image<T>& GetLevel(int i) const { return images[i]; }

        inline float ScaleX(int level) const { return scalesX[level]; }
        inline float ScaleY(int level) const { return scalesY[level]; }

        inline size_t NumLevels() const { return images.size(); }
        
        // The maximum level is the floor(min(log2(w), log2(h)))
        static size_t ComputeMaxLevels(float s, int w, int h)
        {
            if (s >= 1.0f)
            {
                return ~0; // max
            }

            const float invs = 1.0f / s;
            int c = 1;
            float x = 1.0f;
            while (x < w && x < w)
            {
                x *= invs;
                ++c;
            }

            return c;
        }

    private:
        std::vector<Image<T>> images;
        std::vector<float> scalesX;
        std::vector<float> scalesY;
    };

} // namespace sight
