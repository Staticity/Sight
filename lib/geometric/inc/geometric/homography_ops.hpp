#pragma once

#include <Eigen/Core>

#include <linear/vec.hpp>

namespace sight
{
	template <typename S, int N>
	Eigen::Vector<S, N> ToEigen(const Vec<S, N>& v)
	{
		Eigen::Vector<S, N> ev;
		for (int i = 0; i < N; ++i)
		{
			ev(i) = v(i);
		}
		return ev;
	}
	template <typename S, int N>
	Vec<S, N> ToVec(const Eigen::Vector<S, N>& ev)
	{
		Vec<S, N> v;
		for (int i = 0; i < N; ++i)
		{
			v(i) = ev(i);
		}
		return ev;
	}

    template <typename S>
	inline Eigen::Vector3<S> ToHomogeneous(const Eigen::Vector2<S>& v)
	{
		return Eigen::Vector3<S>(v(0), v(1), 1.0);
	}

    template <typename S>
	inline Vec3<S> ToHomogeneous(const Vec2<S>& v)
	{
		return { v(0), v(1), 1.0 };
	}

	template <typename S>
	inline Eigen::Vector2<S> FromHomogeneous(const Eigen::Vector3<S>& v)
	{
		return Eigen::Vector2<S>(v(0) / v(2), v(1) / v(2));
	}

	template <typename S>
	inline Vec2<S> FromHomogeneous(const Vec3<S>& v)
	{
		return { v(0) / v(2), v(1) / v(2) };
	}

	template <typename S>
	inline Eigen::Vector2<S> ApplyHomography(
		const Eigen::Matrix<S, 3, 3>& H,
		const Eigen::Vector2<S>& v)
	{
		return FromHomogeneous<S>(H * ToHomogeneous<S>(v));
	}

	template <typename S>
	inline Vec2<S> ApplyHomography(
		const Eigen::Matrix<S, 3, 3>& H,
		const Vec2<S>& v)
	{
		return FromHomogeneous<S>(H * ToHomogeneous<S>(v));
	}
}