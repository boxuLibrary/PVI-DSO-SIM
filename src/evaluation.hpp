

#pragma once
#include <Eigen/Geometry>
#include <iostream>
#include <vector>

void evaluateTrajectory(std::vector<Eigen::Matrix4d>& Twcs_gt, std::vector<Eigen::Matrix4d>& Twcs_esti);

void evaluatePlane(std::vector<Eigen::Vector4d>& plane_gt, std::vector<Eigen::Vector4d>& plane_est);