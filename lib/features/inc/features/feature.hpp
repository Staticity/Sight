#pragma once

#include <image/image_ops.hpp>

namespace sight
{
    struct Feature
    {
        Feature(
            float x = 0.f,
            float y = 0.f,
            float angle = 0.f,
            float response = 0.f,
            uint8_t level = 0,
            float scaleX = 1.f,
            float scaleY = 1.f)
            : x(x)
            , y(y)
            , level(level)
            , scaleX(scaleX)
            , scaleY(scaleY)
            , angle(angle)
            , response(response)
        {}

        inline float originX() const
        {
            return ScalePx<float>(x, scaleX);
        }

        inline float originY() const
        {
            return ScalePx<float>(y, scaleY);
        }

        bool operator==(const Feature& f) const
        {
            if (response != f.response)
                return false;

            if (level != f.level)
                return false;

            if (scaleX != f.scaleX)
                return false;

            if (scaleY != f.scaleY)
                return false;

            if (y != f.y)
                return false;
            
            if (x != f.x)
                return false;

            if (angle != f.angle)
                return false;

            return true;
        }

        bool operator!=(const Feature& f) const
        {
            return !(*this == f);
        }

        bool operator<(const Feature& f) const
        {
            if (response != f.response)
                return response < f.response;

            if (level != f.level)
                return level < f.level;

            if (scaleX != f.scaleX)
                return scaleX < f.scaleX;

            if (scaleY != f.scaleY)
                return scaleY < f.scaleY;

            if (y != f.y)
                return y < f.y;

            if (x != f.x)
                return x < f.x;
            
            return angle < f.angle;
        }

        bool operator>(const Feature& f) const
        {
            if (response != f.response)
                return response > f.response;

            if (level != f.level)
                return level > f.level;

            if (scaleX != f.scaleX)
                return scaleX > f.scaleX;

            if (scaleY != f.scaleY)
                return scaleY > f.scaleY;

            if (y != f.y)
                return y > f.y;

            if (x != f.x)
                return x > f.x;
            
            return angle > f.angle;
        }

        float x;
        float y;
        float scaleX;
        float scaleY;
        uint8_t level;
        float angle;
        float response;
    };
}