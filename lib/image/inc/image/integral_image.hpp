#pragma once

namespace sight
{
    template <typename S>
    class IntegralImage
    {
    public:

        template <typename T>
        IntegralImage(const Image<T>& im)
            : m_sum(im.w, im.h, 1)
        {
            assert(im.c == 1);


            // Slower implementation.. Could separate it to remove
            // conditionals. (TODO)
            for (int i = 0; i < im.h; ++i)
            {
                auto src = im.row(i);
                auto dst = m_sum.row(i);

                for (int j = 0; j < im.w; ++j)
                {
                    dst[j] = src[j];

                    if (i > 0)
                    {
                        dst[j] += src[j - im.row_step];
                    }

                    if (j > 0)
                    {
                        dst[j] += src[j - 1];
                    }

                    if (i > 0 && j > 0)
                    {
                        dst[j] -= src[j - im.row_step - 1];
                    }
                }
            }
        }

        S Sum(int x, int y, int w, int h)
        {
            --w;
            --h;
            return
                m_sum(x + w, y + h) -
                m_sum(x - 1, y + h) -
                m_sum(x + w, y - 1) +
                m_sum(x - 1, y - 1);
        }

    private:

        Image<T> m_sum;
    };
}
