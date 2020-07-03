#pragma once

#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <opencv2/imgcodecs.hpp>

namespace sight
{
    struct Roi
    {
        Roi(int x, int y, int w, int h)
            : x(x)
            , y(y)
            , w(w)
            , h(h)
        {}

        Roi translate(int dx, int dy) const { return Roi(x + dx, y + dy, w, h); }
        Roi operator+(int v) const { return Roi(x + v, y + v, w, h); }
        Roi operator-(int v) const { return *this + -v; }

        int x;
        int y;
        int w;
        int h;
    };

    template <typename T>
    class Image
    {
        Image(const Image<T>& im, const Roi& roi)
            : w(roi.w)
            , h(roi.h)
            , c(im.c)
            , px(roi.w * roi.h * im.c)
            , row_step(im.row_step)
            , col_step(im.col_step)
            , pixels(im.pixels)
            , sub_image(true)
        {
            pFirst = const_cast<T*>(im.begin()) + (row_step * roi.y) + roi.x * col_step;
            pLast = pFirst + (row_step * roi.h);
        }

    public:
        Image()
            : w(0)
            , h(0)
            , c(0)
            , px(0)
            , row_step(0)
            , col_step(0)
            , sub_image(false)
            , pixels()
        {}

        Image(int width, int height, int channels = 1)
            : w(width)
            , h(height)
            , c(channels)
            , px(width * height * channels)
            , row_step(width * channels)
            , col_step(channels)
            , sub_image(false)
            , pixels(std::make_shared<std::vector<T>>(width * height * channels))
        {
            pFirst = &(*pixels)[0];
            pLast = pFirst + px;
        }

        ~Image()
        {}

        inline T* begin() { return pFirst; }
        inline const T* begin() const { return pFirst; }

        inline T* end() { return pLast; }
        inline const T* end() const { return pLast; }

        inline T* row(int r = 0) { return pFirst + r * row_step; };
        inline const T* row(int r = 0) const { return pFirst + r * row_step; };

        inline T* at(int i, int j) { return pFirst + i * row_step + j * col_step; }
        inline const T* at(int i, int j) const { return pFirst + i * row_step + j * col_step; }

        inline T* at(int i, int j, int ch) { return pFirst + i * row_step + j * col_step + ch; }
        inline const T* at(int i, int j, int ch) const { return pFirst + i * row_step + j * col_step + ch; }

        inline T& operator()(int i) { return begin()[i]; }
        inline const T& operator()(int i) const { return begin()[i]; }

        inline T& operator()(int i, int j) { return *at(i, j); }
        inline const T& operator()(int i, int j) const { return *at(i, j); }

        inline T& operator()(int i, int j, int ch) { return *at(i, j, ch); }
        inline const T& operator()(int i, int j, int ch) const { return *at(i, j, ch); }

        inline Roi GetROI() const { return Roi(0, 0, w, h); }

        Image operator()(const Roi& roi) const
        {
            return Image(*this, roi);
        }

        cv::Mat ToOpenCV() const
        {
            return ToOpenCV(*this);
        }

        static Image Read(const std::string& filepath)
        {
            assert(std::filesystem::exists(filepath));
            return FromOpenCV(cv::imread(filepath));
        }

        static Image FromOpenCV(const cv::Mat& m)
        {
            assert(m.isContinuous());
            Image<T> im(m.cols, m.rows, m.channels());
            im.pixels->assign((T*)m.datastart, (T*)m.dataend);
            return im;
        }

        static cv::Mat ToOpenCV(const Image& im)
        {
            cv::Mat m(im.h, im.w, GetOpenCVType(im));
            
            if (!im.sub_image)
            {
                memcpy(m.data, im.begin(), sizeof(T) * im.px);
            }
            else
            {
                auto dst = m.ptr<T>();
                for (int i = 0; i < im.h; ++i)
                {
                    auto srcRow = im.row(i);
                    for (int j = 0; j < im.w; ++j)
                    {
                        for (int k = 0; k < im.c; ++k)
                        {
                            *(dst++) = *(srcRow++);
                        }
                    }
                }
            }
            
            return m;
        }

        int w;
        int h;
        int c;
        int row_step;
        int col_step;
        int px;

        T* pFirst;
        T* pLast;

    private:

        static int GetOpenCVType(const Image<T>& im);

        bool sub_image;
        std::shared_ptr<std::vector<T>> pixels;
    };

    template <typename T>
    std::ostream& operator<<(std::ostream& o, const Image<T>& im)
    {
        for (int ch = 0; ch < im.c; ++ch)
        {
            if (im.c > 1)
                o << "Channel: " << ch << '\n';
            for (int i = 0; i < im.h; ++i)
            {
                for (int j = 0; j < im.w; ++j)
                {
                    o << std::setw(5) << +im(i, j, ch) << ' ';
                }
                o << '\n';
            }
        }

        return o;
    }

} // namespace sight