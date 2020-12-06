#pragma once 

#include <vector>
#include <string>

namespace sight
{
	struct VideoSourceInfo
	{
		bool valid{ false };
		size_t index{ 0 };
		std::string name;
		std::string description;
		std::string identifier;
	};

	VideoSourceInfo GetVideoSource(int connectedCameraIdx);
	VideoSourceInfo GetVideoSource(const std::string& identifier);

	std::vector<VideoSourceInfo> EnumerateVideoSources();
}
