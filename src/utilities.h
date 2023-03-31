//
// Created by hyj on 18-1-19.
//

#ifndef IMUSIMWITHPOINTLINE_UTILITIES_H
#define IMUSIMWITHPOINTLINE_UTILITIES_H

#include "imu.h"
#include "random.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <vector>
#include <fstream>

// save 3d points to file
void save_points(std::string filename, std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > points);
void LoadPoints(std::string filename, std::vector<Eigen::Vector3d>& pts);

// save 3d points and it's obs in image
void save_features(std::string filename,
                   std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > points,
                   std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > features);

// save line obs
void save_lines(std::string filename,
                std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > features);

void LoadPointObs(std::string filename, std::vector<Eigen::Vector3d>& pts ,std::vector<Eigen::Vector3d>& obs);
void LoadPose(std::string filename, std::vector<MotionData>& pose);

// save imu body data
void save_Pose(std::string filename, std::vector<MotionData> pose);

// save pose as TUM style
void save_Pose_asTUM(std::string filename, std::vector<MotionData> pose);


// math tools
template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymmetric(const Eigen::MatrixBase<Derived> &q)
{
    Eigen::Matrix<typename Derived::Scalar, 3, 3> ans;
    ans << typename Derived::Scalar(0), -q(2), q(1),
            q(2), typename Derived::Scalar(0), -q(0),
            -q(1), q(0), typename Derived::Scalar(0);
    return ans;
}
template<class bidiiter> //Fisher-Yates shuffle
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random)
{
    size_t left = std::distance(begin, end);

    while (num_random--)
    {
        bidiiter r = begin;
        std::advance(r, xrand() % left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }

    return begin;
}

#endif //IMUSIMWITHPOINTLINE_UTILITIES_H


