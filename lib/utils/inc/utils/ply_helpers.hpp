#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>

#include <linear/vec.hpp>

namespace sight
{
    template <typename S>
    void WritePlyPointCloud(
        const std::string& filepath,
        const std::vector<Vec3<S>>& points,
        const std::vector<Vec3<uint8_t>>& colors = std::vector<Vec3<uint8_t>>())
    {
        std::ofstream out(filepath);

        out << "ply\n";
        out << "format ascii 1.0\n";
        out << "element vertex " << points.size() << "\n";
        out << "property float x\n";
        out << "property float y\n";
        out << "property float z\n";
        out << "property uchar red\n";
        out << "property uchar green\n";
        out << "property uchar blue\n";
        out << "end_header\n";
        for (int i = 0; i < points.size(); ++i)
        {
            out << points[i][0] << ' ';
            out << points[i][1] << ' ';
            out << points[i][2] << ' ';
            if (i < colors.size())
            {
                out << colors[i][0] << ' ';
                out << colors[i][1] << ' ';
                out << colors[i][2] << '\n';
            }
            else
            {
                out << "128 128 128\n";
            }
        }
    }
}