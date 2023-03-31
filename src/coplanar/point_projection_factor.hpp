#pragma once
#include <ceres/ceres.h>
#include <sophus/se3.hpp>
#include <Eigen/Dense>
#include <ceres/rotation.h>
#include "utility.hpp"
#include "tic_toc.h"

class ProjectionPointInverseDepthFactor : public ceres::SizedCostFunction<2, 7, 7, 1>
{
public:
    ProjectionPointInverseDepthFactor(const Eigen::Vector3d& ob0, const Eigen::Vector3d& obi);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    static Eigen::Matrix2d sqrt_info;
    static double sum_t;
    static Eigen::Quaterniond qic;
    static Eigen::Vector3d tic;

    Eigen::Vector3d point_obi_;
    Eigen::Vector3d point_ob0_;

};


class PointOnPlaneFactor : public ceres::SizedCostFunction<1, 4, 7, 1>
{
public:
    PointOnPlaneFactor(const Eigen::Vector3d &_pts_i);

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);
    static double sqrt_info;
    static Eigen::Vector3d tic;
    static Eigen::Quaterniond qic;

    Eigen::Vector3d pts_i;


};

class PointOnPlaneProjectionFactor : public ceres::SizedCostFunction<2, 4, 7, 7>
{
public:
    PointOnPlaneProjectionFactor(
            const Eigen::Vector3d&_ob0,
            const Eigen::Vector3d& _obi);

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);
    static Eigen::Matrix2d sqrt_info;
    static Eigen::Quaterniond qic;
    static Eigen::Vector3d tic;
    Eigen::Vector3d ob0;
    Eigen::Vector3d obi;
};




class PointOnPlaneProjectionPriorFactorWall : public ceres::SizedCostFunction<3, 4>
{
public:
    PointOnPlaneProjectionPriorFactorWall(
            const Eigen::Vector4d& plane_prior_);

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);
    static Eigen::Matrix3d sqrt_info;
    Eigen::Vector4d plane_prior;

};


class PointOnPlaneProjectionPriorFactorGround: public ceres::SizedCostFunction<3, 4>
{
public:
    PointOnPlaneProjectionPriorFactorGround(
            const Eigen::Vector4d& plane_prior_);

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);
    static Eigen::Matrix3d sqrt_info;
    Eigen::Vector4d plane_prior;

};