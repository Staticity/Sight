#include <iostream>
#include <numeric>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>

#include <sophus/se3.hpp>

#include <opencv2/opencv.hpp>

#include <image/image.hpp>
#include <image/image_gen.hpp>
#include <image/image_ops.hpp>
#include <image/pyramid.hpp>
#include <features/fast.hpp>
#include <features/orb.hpp>
#include <geometric/homography_estimation.hpp>

using std::cout;
using std::endl;
using std::cin;

using namespace cv;

void MatchPair(const std::string& filepath1, const std::string& filepath2)
{
    using namespace sight;
    const auto im1 = Image<uint8_t>::Read(filepath1);
    const auto im2 = Image<uint8_t>::Read(filepath2);

    const auto pyr1 = Pyramid<uint8_t>(ToGrayscale(im1), 10, .8f);
    const auto pyr2 = Pyramid<uint8_t>(ToGrayscale(im2), 10, .8f);
    
    // for (int i = 0; i < pyr1.NumLevels(); ++i)
    // {
    //     imshow("Pyramid Level: " + std::to_string(i), pyr1.GetLevel(i).ToOpenCV());
    //     waitKey(0);
    // }

    const int nMaxFeatures = 2000;
    const auto feats1 = sight::FindFAST(pyr1);
    const auto descs1 = sight::ComputeORB(pyr1, feats1, nMaxFeatures);

    const auto feats2 = sight::FindFAST(pyr2);
    const auto descs2 = sight::ComputeORB(pyr2, feats2, nMaxFeatures);

    // Brute force matching
    std::vector<int> index2FromIndex1(descs1.size());
    std::vector<int> bestDistance(descs1.size(), 256);
    std::vector<float> bestRatio(descs1.size(), 1.0f);
    for (int i = 0; i < descs1.size(); ++i)
    {
        if (descs1[i].feat.level == 0) continue;
        const auto& d1 = descs1[i];
        for (int j = 0; j < descs2.size(); ++j)
        {
            const auto& d2 = descs2[j];
            const auto dist = d1.distance(d2);
            if (dist < bestDistance[i])
            {
                bestRatio[i] = dist / float(bestDistance[i]);
                bestDistance[i] = dist;
                index2FromIndex1[i] = j;
            }
        }
    }

    // Fit a homography via RANSAC with the 'good' matches
    const float goodRatio = .75f;
    EigenList<Eigen::Vector2d> pts1, pts2;
    pts1.reserve(descs1.size());
    pts2.reserve(descs2.size());
    for (int i = 0; i < descs1.size(); ++i)
    {
        const int j = index2FromIndex1[i];
        const auto& f1 = descs1[i].feat;
        const auto& f2 = descs2[j].feat;

        if (bestRatio[i] < goodRatio)
        {
            pts1.push_back(Eigen::Vector2d(f1.originX(), f1.originY()));
            pts2.push_back(Eigen::Vector2d(f2.originX(), f2.originY()));
        }
    }

    std::vector<int> indices;
    // std::vector<int> indices(pts1.size());
    // std::iota(indices.begin(), indices.end(), 0);
    RansacHomography(pts1, pts2, 200, 15.0, indices);
    
    auto merge = sight::MergeHorizontal(im1, im2).ToOpenCV();

    RNG rng(12345);
    for (int i : indices)
    {
        const Scalar color(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));

        const auto p1 = cv::Point2d(pts1[i].x(), pts1[i].y());
        const auto p2 = cv::Point2d(im1.w + pts2[i].x(), pts2[i].y());
        cv::line(merge, p1, p2, color);
    }

    cv::resize(merge, merge, Size(), .75, .75);
    imshow("Matching", merge);
    waitKey(0);
    exit(0);
}

int main(int argc, char** argv)
{
    MatchPair(
        "C:\\Users\\Jaime\\source\\repos\\Sight\\assets\\homography\\aerial\\0.jpg",
        "C:\\Users\\Jaime\\source\\repos\\Sight\\assets\\homography\\aerial\\1.jpg"
    );


    // auto im = sight::CreateCorner<uint8_t>(499, M_PI_2, 0, 255);

    const std::string filepath =
        "C:\\Users\\Jaime\\Downloads\\simple.jpg";

    auto image = sight::Image<uint8_t>::Read(filepath);
    auto gray = sight::ToGrayscale(image);

    const sight::Pyramid<uint8_t> pyr(gray, 5, .8f);
    auto feats = sight::FindFAST(pyr);
    auto descs = sight::ComputeORB(pyr, feats, 500);

    std::vector<cv::Mat> levelDrawings(pyr.NumLevels());
    for (int i = 0; i < pyr.NumLevels(); ++i)
    {
        levelDrawings[i] = pyr.GetLevel(i).ToOpenCV();
    }

#if 0
    for (int i = 0; i < pyr.NumLevels(); ++i)
    {
        auto temp = pyr.GetLevel(i).ToOpenCV();
        for (const auto& feat : feats)
        {
            if (feat.level == i)
            {
                cv::circle(temp, cv::Point2d(feat.x, feat.y), 3, cv::Scalar(0, 0, 255));
            }
        }
        imshow("level: " + std::to_string(i), temp);
    }
#endif    

    cv::Mat mat = image.ToOpenCV();
    for (const auto& desc : descs)
    {
        const auto& feat = desc.feat;
        const auto center = cv::Point2d(feat.originX(), feat.originY());
        const int radius = 4;
        const int size = int(round(3 * feat.scaleX + 1));
        const auto color = cv::Scalar(0, 0, 255);

        const cv::Point2d dir(cos(feat.angle), sin(feat.angle));

        cv::circle(mat, center, size, color);
        cv::arrowedLine(mat, center, center + dir * size, color);

        cv::circle(levelDrawings[feat.level], Point2d(feat.x, feat.y), 3, color);
        cv::arrowedLine(levelDrawings[feat.level], Point2d(feat.x, feat.y), Point2d(feat.x, feat.y) + dir * 3, color);
    }

    int i = 0;
    for (auto& drawing : levelDrawings)
    {
        imshow("Level Drawing: " + std::to_string(i++), drawing);
    }

    imwrite("debug_image.png", mat);
    imshow("test", mat);
    waitKey(0);

    return 0;
}