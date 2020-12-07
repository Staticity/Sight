#pragma once

#include <cmath>
#include <functional>
#include <algorithm>
#include <array>

#include <linear/vec_ops.hpp>

namespace sight
{
    template <typename S, int N>
    class Vec;

    template <typename S>
    using Vec2 = Vec<S, 2>;

    template <typename S>
    using Vec3 = Vec<S, 3>;

    template <typename S>
    using Vec4 = Vec<S, 4>;

    template <typename S, int N>
    class Vec
    {
    public:
        Vec()
            : v()
        {}

        Vec(S x)
            : v()
        {
            Fill(x);
        }

        Vec(std::initializer_list<S> l)
        {
            std::copy(l.begin(), l.end(), v.begin());
        }

        inline S& operator[](int i) { return v[i]; }
        inline const S& operator[](int i) const { return v[i]; }

        inline S& operator()(int i) { return v[i]; }
        inline const S& operator()(int i) const { return v[i]; }

        Vec operator-() const
        {
            Vec neg;
            for (int i = 0; i < N; ++i)
            {
                neg[i] = -v[i];
            }
            return neg;
        }

        Vec operator*(const S s) const
        {
            Vec res;
            for (int i = 0; i < N; ++i)
            {
                res[i] = s * v[i];
            }
            return res;
        }

        Vec& operator*=(const S s)
        {
            for (int i = 0; i < N; ++i)
            {
                v[i] *= s;
            }
            return *this;
        }
        
        Vec operator+(const Vec& u) const
        {
            return Apply(u, std::plus<S>());
        }

        Vec& operator+=(const Vec& u)
        {
            *this = *this + u;
            return *this;
        }

        Vec operator-(const Vec& u) const
        {
            return Apply(u, std::minus<S>());
        }

        Vec& operator-=(const Vec& u)
        {
            *this = *this - u;
            return *this;
        }

        Vec Normalized() const
        {
            Vec x = *this;
            x.Normalize();
            return x;
        }

        Vec3<S> Cross(const Vec3<S>& u) const
        {
            Vec3<S> vxu;
            Cross3(&v[0], &u[0], &vxu[0]);
            return vxu;
        }

        void Normalize()
        {
            const S invnorm = S(1) / Norm();
            for (int i = 0; i < N; ++i)
            {
                v[i] *= invnorm;
            }
        }

        S Norm() const
        {
            return sqrt(SquaredNorm());
        }

        S SquaredNorm() const
        {
            return Dot(*this);
        }

        S Dot(const Vec& u) const
        {
            S x = S();
            for (int i = 0; i < N; ++i)
            {
                x += v[i] * u[i];
            }
            return x;
        }

        void Zero()
        {
            return Fill(S());
        }

        void Fill(S x)
        {
            for (int i = 0; i < N; ++i)
            {
                v[i] = x;
            }
        }

    private:
        
        template <typename Op>
        inline Vec Apply(const Vec& u, const Op& op) const
        {
            Vec x;
            for (int i = 0; i < N; ++i)
            {
                x[i] = op(v[i], u[i]);
            }
            return x;
        }

        std::array<S, N> v;
    };

    template <typename S, int N>
    inline Vec<S, N> operator*(S s, const Vec<S, N>& v)
    {
        return v * s;
    }
}