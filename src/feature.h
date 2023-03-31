//
// Created by xubo on 22-11-24.
//

#ifndef VIO_INIT_SYS_FEATURE_H
#define VIO_INIT_SYS_FEATURE_H
#include <fstream>
#include <random>
#include <unordered_map>
#include <map>
#include <Eigen/Core>
using Point = Eigen::Vector4d;
using Points = std::vector<Point, Eigen::aligned_allocator<Point> >;
using FeatureID = int;
using PlaneID = int;
using TimeFrameId = int;

class FeaturePerFrame {
public:
    FeaturePerFrame() = default;

    FeaturePerFrame(const Eigen::Vector3d &normalPt, const Eigen::Vector2d &uvPt) {
        normalpoint = normalPt;
        uv = uvPt;
    }

    Eigen::Vector3d normalpoint;
    Eigen::Vector2d uv;
    double z;
    bool is_used;
    double parallax;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    double dep_gradient;
};

class SFMFeature {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SFMFeature() = default;

    SFMFeature(FeatureID _feature_id, PlaneID _plane_id, TimeFrameId _start_frame)
            : feature_id(_feature_id), plane_id(_plane_id), kf_id(_start_frame) {}

    FeatureID feature_id{};
    TimeFrameId kf_id{};
    PlaneID plane_id{};

    bool state = false; // sfm使用：被三角化: true, 未被三角化：false
    Eigen::Vector3d p3d; // sfm使用：世界坐标系下的3D位置
    double Idepth;
    double Idepth_opt;
    double Idepth_gt;
    // TODO: 观测值
    std::map<TimeFrameId, FeaturePerFrame> obs;
};
#endif //VIO_INIT_SYS_FEATURE_H
