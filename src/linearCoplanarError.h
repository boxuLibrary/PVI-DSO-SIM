//
// Created by xubo on 22-11-22.
//

#ifndef VIO_INIT_SYS_LINEARCOPLANARERROR_H
#define VIO_INIT_SYS_LINEARCOPLANARERROR_H

#include <sophus/se3.hpp>
#include "param.h"
#include "utilities.h"
#include "feature.h"

class BundleAdjustmentBase {
public:
    BundleAdjustmentBase(const Param &params,
                         std::unordered_map<FeatureID, SFMFeature> &points_obs_,
                         std::vector<Eigen::Matrix4d> &Twcs_,
                         std::vector<Eigen::Vector4d> &planes_
    ) : width(params.image_w),
        height(params.image_h),
        Ric(params.R_bc),
        tic(params.t_bc),
        optimization_num(params.optimization_num),
        points_obs(points_obs_),
        Twcs(Twcs_),
        planes(planes_) {


        Eigen::Matrix3d K_;
        K_.setIdentity();
        K_(0, 0) = params.fx;
        K_(0, 2) = params.cx;
        K_(1, 1) = params.fy;
        K_(1, 2) = params.cy;
        K = K_;

        constraint_type = params.constraint_type;

        Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
        Tic.block(0, 0, 3, 3) = Ric;
        Tic.block(0, 3, 3, 1) = tic;
        for (int i = 0; i < Twcs.size(); i++) {
            Twis.emplace_back(Twcs[i] * Tic.inverse());
        }

        planes_prior = planes;

    }

    ~BundleAdjustmentBase() = default;

    void triangulate();

    void add_pose_noise();

    void add_wall_plane_noise();

    void loosely_coupled_method();

    void tightly_coupled_method();

    void point_method();

    double updateCoplanarPoint(const Eigen::Vector3d &lineStart, const Eigen::Vector3d &lineEnd,
                               const Eigen::Vector4d &plane_para_host);


public:
    Eigen::Matrix3d K;
    int width;
    int height;
    int constraint_type;
    int optimization_num;
    std::unordered_map<FeatureID, SFMFeature> points_obs;
    std::vector<Eigen::Matrix4d> Twcs;
    std::vector<Eigen::Matrix4d> Twcs_opt;
    std::vector<Eigen::Matrix4d> Twis;
    std::vector<Eigen::Matrix4d> Twis_opt;
    std::vector<Eigen::Vector4d> planes;
    std::vector<Eigen::Vector4d> planes_prior;
    std::vector<Eigen::Vector4d> planes_opt;
    Eigen::Matrix3d Ric;
    Eigen::Vector3d tic;



};

#endif //VIO_INIT_SYS_LINEARCOPLANARERROR_H
