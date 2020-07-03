#pragma once

#include <Eigen/Eigen>
#include <vector>
#include <map>

namespace sight
{

template <typename T>
using EigenList = std::vector<T, Eigen::aligned_allocator<T>>;

template <typename K, typename V>
using EigenMap = std::map<
	K, V,
	std::less<std::pair<K, V>>,
	Eigen::aligned_allocator<std::pair<K, V>>>;

} // namespace sight