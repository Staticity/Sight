#pragma once

#include <numeric>
#include <random>

#include <geometric/homography_ops.hpp>
#include <Eigen/SVD>

#include <utils/eigentypes.hpp>

namespace sight
{

	template <typename S>
	inline Eigen::Matrix3<S> InvertScaleAndTrans(const Eigen::Matrix3<S>& T)
	{
		// T is of the form:
		//
		// [s 0 sx]
		// [0 s sy]
		// [0 0  1]
		//
		// Which has a closed form of the following:
		//
		// [1/s  0  -x]
		// [ 0  1/s -y]
		// [ 0   0   1]
		const auto s = (T(0, 0) + T(1, 1)) / 2.0; // average both `s` in case..
		const auto si = 1.0 / s;
		const auto x = T(0, 2) * si;
		const auto y = T(1, 2) * si;

		Eigen::Matrix3<S> Tinv;
		Tinv <<
			si,  0, -x,
			0, si, -y,
			0,  0,  1;
		return Tinv;
	}

	template <typename S>
	Eigen::Matrix<S, 3, 3> NormalizationMatrix(
		const EigenList<Eigen::Vector2<S>>& vs)
	{
		// Compute mean
		const size_t n = vs.size();

		Eigen::Vector2<S> mean(0, 0);
		for (const auto& v : vs)
			mean += v;
		mean *= (1.0 / n);

		// Compute average distance from mean
		S mean_mag = 0.0;
		for (const auto& v : vs)
			mean_mag += (v - mean).norm();
		mean_mag /= n;
		
		// We would like the average distance to be
		// the sqrt(2), meaning average of 1 in each
		// component.
		S s = sqrt(2.0) / mean_mag;

		// Compute the normalization matrix:
		//
		// [s 0 0] [1 0 -tx] [x]
		// [0 s 0] [0 1 -ty] [y]
		// [0 0 1] [0 0   1] [1]
		//
		// Which combines into:
		//
		// [s 0 -stx]
		// [0 s -sty]
		// [0 0    1]
		
		Eigen::Matrix<S, 3, 3> T;
		T <<
			s, 0, -s * mean(0),
			0, s, -s * mean(1),
			0, 0, 1;
		return T;
	}

	template <typename S>
	Eigen::Matrix<S, 3, 3> FitHomography(
		const EigenList<Eigen::Vector2<S>>& src,
		const EigenList<Eigen::Vector2<S>>& dst,
		bool normalize=true)
	{
		if (src.size() != dst.size())
			return Eigen::Matrix<S, 3, 3>::Zero();

		// (H * src{i}) x dst{i} = 0
		// [a b c] [x]   [u]
		// [d e f] [y] = [v]
		// [g h i] [z] = [w]
		//
		// set src{i} as [x, y, z]
		// set dst{i} as [u, v, w]
		//
		// Vectorize the homography to solve:
		const size_t n = src.size();

		Eigen::Matrix3<S> T_src; 
		Eigen::Matrix3<S> T_dst;

		if (normalize)
		{
			T_src = NormalizationMatrix<S>(src);
			T_dst = NormalizationMatrix<S>(dst);
		}

		Eigen::Matrix<S, -1, 9> A;
		A.resize(2 * n, 9);

		for (int i = 0; i < n; ++i)
		{
			S x, y, z;
			S u, v, w;

			// Set convenience variables depending
			// on whether we normalize or not.
			if (normalize)
			{
				const Eigen::Vector3d Ts = T_src * ToHomogeneous<S>(src[i]);
				const Eigen::Vector3d Td = T_dst * ToHomogeneous<S>(dst[i]);

				x = Ts(0);
				y = Ts(1);
				z = Td(2);

				u = Td(0);
				v = Td(1);
				w = Td(2);
			}
			else
			{
				x = src[i](0);
				y = src[i](1);
				z = S(1);

				u = dst[i](0);
				v = dst[i](1);
				w = S(1);
			}
			
			// Duplicated work
			const auto xw = x * w;
			const auto yw = y * w;
			const auto zw = z * w;

			// Set rows
			A.row(i * 2 + 0) << 0, 0, 0, xw, yw, zw, -x*v, -y*v, -z*v;
			A.row(i * 2 + 1) << -xw, -yw, -zw, 0, 0, 0, x*u, y*u, z*u;
		}

		Eigen::JacobiSVD<Eigen::Matrix<S, -1, 9>> svd(A, Eigen::ComputeFullV);
		const auto x = svd.matrixV().col(8);
		Eigen::Matrix<S, 3, 3> H;
		H <<
			x(0), x(1), x(2),
			x(3), x(4), x(5),
			x(6), x(7), x(8);

		if (normalize)
		{
			const auto T_dst_inv = InvertScaleAndTrans<S>(T_dst);
			H = T_dst_inv * H * T_src;
		}

		return H;
	}

	template <typename S>
	Eigen::Matrix3<S> RansacHomography(
		const EigenList<Eigen::Vector2<S>>& src,
		const EigenList<Eigen::Vector2<S>>& dst,
		int nIterations,
		S inlierThreshold,
		std::vector<int>& inlierIndices)
	{
		if (src.size() != dst.size() || src.size() < 4)
		{
			return Eigen::Matrix3<S>::Zero();
		}

		const S inlierThresholdSq = pow(inlierThreshold, 2);
		Eigen::Matrix3<S> H;

		const int N = 4;
		std::vector<int> sample(N);
		std::vector<int> indices(src.size());
		std::iota(indices.begin(), indices.end(), 0);

		for (int it = 0; it < nIterations; ++it)
		{
			// Select 4 random indices
			std::sample(
				indices.begin(),
				indices.end(),
				sample.begin(),
				N,
				std::mt19937{std::random_device{}()});
			
			// Store the sample we've selected
			EigenList<Eigen::Vector2<S>> sample_src(N);
			EigenList<Eigen::Vector2<S>> sample_dst(N);
			for (int i = 0; i < N; ++i)
			{
				sample_src[i] = src[sample[i]];
				sample_dst[i] = dst[sample[i]];
			}

			// Count the number of inliers this homography creates
			const Eigen::Matrix3<S> H_guess = FitHomography(sample_src, sample_dst, true);
			std::vector<int> inlierIndices_guess;
			for (int i = 0; i < src.size(); ++i)
			{
				const auto& actual = dst[i];
				const auto pred = ApplyHomography(H_guess, src[i]);
				const auto errSq = (actual - pred).squaredNorm();

				if (errSq <= inlierThresholdSq)
				{
					inlierIndices_guess.push_back(i);
				}
			}

			// Ignore homographies which can't map the
			// set of 4 correspondences.
			if (inlierIndices_guess.size() < N) continue;

			// This guess had more inilers, let's remember it.
			if (inlierIndices_guess.size() > inlierIndices.size())
			{
				inlierIndices = inlierIndices_guess;
				H = H_guess;
			}
		}

		return H;
	}

} // namespace sight
