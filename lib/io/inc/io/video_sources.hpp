#pragma once 

#include <vector>
#include <string>

namespace sight
{
	struct VideoSourceInfo
	{
		size_t index{ 0 };
		std::string name;
		std::string description;
		std::string identifier;
	};

	std::vector<VideoSourceInfo> EnumerateVideoSources();
}
