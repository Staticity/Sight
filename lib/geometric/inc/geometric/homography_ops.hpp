#pragma once

#include <Eigen/Core>

namespace sight
{
    template <typename S>
	inline Eigen::Vector3<S> ToHomogeneous(const Eigen::Vector2<S>& v)
	{
		return Eigen::Vector3<S>(v(0), v(1), 1.0);
	}

	template <typename S>
	inline Eigen::Vector2<S> FromHomogeneous(const Eigen::Vector3<S>& v)
	{
		return Eigen::Vector2<S>(v(0) / v(2), v(1) / v(2));
	}

	template <typename S>
	inline Eigen::Vector2<S> ApplyHomography(
		const Eigen::Matrix<S, 3, 3>& H,
		const Eigen::Vector2<S>& v)
	{
		return FromHomogeneous<S>(H * ToHomogeneous<S>(v));
	}

}