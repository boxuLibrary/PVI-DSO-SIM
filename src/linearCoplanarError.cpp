//
// Created by xubo on 22-11-22.
//
#include "linearCoplanarError.h"
#include "coplanar/pose_local_parameterization.hpp"
#include "coplanar/plane_local_parameterization.hpp"
#include "coplanar/point_projection_factor.hpp"

void BundleAdjustmentBase::triangulate() {

    for (auto &pt: points_obs) {
        Eigen::MatrixXd A(pt.second.obs.size() * 2, 4);
        TimeFrameId frame_id0 = pt.second.kf_id;
        Eigen::Matrix4d Twc0 = Twcs[frame_id0];
        int index = 0;
        for (const auto &frame: pt.second.obs) {

            Eigen::Matrix4d P = Twcs[frame.first].inverse() * Twc0;
            Eigen::Vector3d f = frame.second.normalpoint;
            A.row(index++) = f[0] * P.row(2) - f[2] * P.row(0);
            A.row(index++) = f[1] * P.row(2) - f[2] * P.row(1);
        }

        assert(index == A.rows());
        Eigen::Vector4d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(A,
                                                                  Eigen::ComputeThinV).matrixV().rightCols<1>();
        double svd_method = svd_V[2] / svd_V[3];
        pt.second.state = true;
        pt.second.Idepth = 1 / svd_method;
    }

}


void BundleAdjustmentBase::add_pose_noise() {

    double max_nt = 0.3;
    double max_nq = 5. * M_PI / 180.;

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> nt(0., max_nt);
    std::normal_distribution<double> nq(0., max_nq);

    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.block(0, 3, 3, 1) = tic;
    Tic.block(0, 0, 3, 3) = Ric;


    for (int index = 0; index < Twcs.size(); ++index) {

        if (index < 10) continue;

        auto &Twc = Twcs[index];
        auto &Twi = Twis[index];
        Eigen::Matrix4d Tcc = Eigen::Matrix4d::Identity();
        Eigen::Vector3d tcc = Eigen::Vector3d(nt(generator), nt(generator), nt(generator));
        Eigen::Matrix3d Rcc = Eigen::Matrix3d::Identity();
        Rcc = Eigen::AngleAxisd(nq(generator), Eigen::Vector3d::UnitZ())
              * Eigen::AngleAxisd(nq(generator), Eigen::Vector3d::UnitY())
              * Eigen::AngleAxisd(nq(generator), Eigen::Vector3d::UnitX());

        Tcc.block(0, 0, 3, 3) = Rcc;
        Tcc.block(0, 3, 3, 1) = tcc;
        Twc = Twc * Tcc;
        Twi = Twc * Tic.inverse();
    }
}

void BundleAdjustmentBase::add_wall_plane_noise() {

    double max_distance = 0.3;
    double max_normal = 5 * M_PI / 180.;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> n_distance(0., max_distance);
    std::normal_distribution<double> n_normal(0., max_normal);

    for (auto &plane: planes) {

        Eigen::Vector3d norm_w = plane.head<3>();

        double theta = std::atan2(norm_w.y(), norm_w.x());
        double phi = std::asin(norm_w.z());

        theta += n_normal(generator);
        phi += n_normal(generator);

        Eigen::Vector3d n(cos(phi) * cos(theta), cos(phi) * sin(theta), sin(phi));

        double d = plane.w();
        d += n_distance(generator);

        plane.head<3>() = n;
        plane(3) = d;
    }
}


void BundleAdjustmentBase::point_method() {

    // weight of point constraints
    ProjectionPointInverseDepthFactor::sqrt_info = 460.0 / 1.5 * Eigen::Matrix2d::Identity();
    ProjectionPointInverseDepthFactor::qic = Eigen::Quaterniond(Ric);
    ProjectionPointInverseDepthFactor::tic = tic;

    ceres::Problem problem;

    // points parameters
    for (auto &pt: points_obs) {
        pt.second.Idepth_opt = pt.second.Idepth;
        if (pt.second.state) {
            problem.AddParameterBlock(&pt.second.Idepth_opt, 1);
        }
    }

    // pose parameters
    std::vector<Eigen::Matrix<double, 7, 1>> para_Pose(Twis.size());
    for (int i = 0; i < Twis.size(); ++i) {
        Eigen::Matrix3d Rot = Twis[i].block(0, 0, 3, 3);
        Eigen::Vector3d Tran = Twis[i].block(0, 3, 3, 1);
        Eigen::Quaterniond qua(Rot);

        para_Pose[i][0] = Tran(0);
        para_Pose[i][1] = Tran(1);
        para_Pose[i][2] = Tran(2);
        para_Pose[i][3] = qua.x();
        para_Pose[i][4] = qua.y();
        para_Pose[i][5] = qua.z();
        para_Pose[i][6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i].data(), 7, local_parameterization);

        // fix 10 frames
        if (i < 10)
            problem.SetParameterBlockConstant(para_Pose[i].data());
    }

    for (auto &pt: points_obs) {

        TimeFrameId first_frame_id = pt.second.kf_id;

        auto &obs = pt.second.obs;

        // add point constraints
        {
            for (const auto &ob: obs) {
                if (ob.first == first_frame_id) continue;

                ceres::LossFunction *loss_function = nullptr;
                loss_function = new ceres::CauchyLoss(1.);
                auto cost_function = new ProjectionPointInverseDepthFactor(obs.at(first_frame_id).normalpoint,
                                                                           ob.second.normalpoint);

                problem.AddResidualBlock(cost_function, loss_function, para_Pose[first_frame_id].data(),
                                         para_Pose[ob.first].data(), &pt.second.Idepth_opt);

            }
        }
    }


    ceres::Solver::Options options;
    // options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;  // LEVENBERG_MARQUARDT  DOGLEG
    // options.linear_solver_type = ceres::SPARSE_SCHUR; // SPARSE_NORMAL_CHOLESKY  or SPARSE_SCHUR
    options.max_num_iterations = optimization_num;
    options.minimizer_progress_to_stdout = false;

    TicToc solver_time;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // std::cout << "loose solver time: " <<  solver_time.toc() << std::endl;
    // info.opti_step = summary.iterations.size();
    std::cout << summary.BriefReport() << std::endl;

    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.block(0, 3, 3, 1) = tic;
    Tic.block(0, 0, 3, 3) = Ric;

    Twis_opt.resize(Twcs.size());
    Twcs_opt.resize(Twcs.size());

    for (int i = 0; i < Twcs.size(); ++i) {
        Eigen::Quaterniond qua(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]);
        Eigen::Matrix3d rot = qua.toRotationMatrix();
        Eigen::Vector3d tran(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);

        Eigen::Matrix4d Twi = Eigen::Matrix4d::Identity();
        Twi.block(0, 3, 3, 1) = tran;
        Twi.block(0, 0, 3, 3) = rot;
        Twcs_opt[i] = Twi * Tic;
        Twis_opt[i] = Twi;
    }
}

void BundleAdjustmentBase::loosely_coupled_method() {

    // weight of two constraints
    ProjectionPointInverseDepthFactor::sqrt_info = 460.0 / 1.5 * Eigen::Matrix2d::Identity();
    PointOnPlaneFactor::sqrt_info = 460.0 / 0.5;
    PointOnPlaneProjectionPriorFactorWall::sqrt_info = 1e5 * Eigen::Matrix3d::Identity();

    PointOnPlaneProjectionPriorFactorGround::sqrt_info = 1e7 * Eigen::Matrix3d::Identity();

    ProjectionPointInverseDepthFactor::qic = Eigen::Quaterniond(Ric);
    ProjectionPointInverseDepthFactor::tic = tic;

    PointOnPlaneFactor::qic = Eigen::Quaterniond(Ric);
    PointOnPlaneFactor::tic = tic;

    ceres::Problem problem;

    // points parameters
    for (auto &pt: points_obs) {

        pt.second.Idepth_opt = pt.second.Idepth;

        if (pt.second.state) {
            problem.AddParameterBlock(&pt.second.Idepth_opt, 1);
        }
    }

    planes_opt.clear();
    planes_opt = planes;

    // plane parameters
    for (int i = 0; i < planes_opt.size(); i++) {
        ceres::LocalParameterization *local_parameterization;
        local_parameterization = new mesh::WallLocalParameterization();
        problem.AddParameterBlock(planes_opt[i].data(), 4, local_parameterization);

    }

    // pose parameters
    std::vector<Eigen::Matrix<double, 7, 1>> para_Pose(Twis.size());
    for (int i = 0; i < Twis.size(); ++i) {
        Eigen::Matrix3d Rot = Twis[i].block(0, 0, 3, 3);
        Eigen::Vector3d Tran = Twis[i].block(0, 3, 3, 1);
        Eigen::Quaterniond qua(Rot);

        para_Pose[i][0] = Tran(0);
        para_Pose[i][1] = Tran(1);
        para_Pose[i][2] = Tran(2);
        para_Pose[i][3] = qua.x();
        para_Pose[i][4] = qua.y();
        para_Pose[i][5] = qua.z();
        para_Pose[i][6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i].data(), 7, local_parameterization);

        // fix 10 frames
        if (i < 10)
            problem.SetParameterBlockConstant(para_Pose[i].data());
    }

    for (auto &pt: points_obs) {

        // 根据特征点的点号找到落在的平面位置
        auto plane_id = pt.second.plane_id;
        // 第一帧对应的帧号
        TimeFrameId first_frame_id = pt.second.kf_id;

        auto &obs = pt.second.obs;

        // add plane constraints
        {
            ceres::LossFunction *loss = nullptr;
            loss = new ceres::CauchyLoss(1.);
            auto cost_plane = new PointOnPlaneFactor(obs.at(first_frame_id).normalpoint);
            problem.AddResidualBlock(cost_plane, loss, planes_opt[plane_id].data(), para_Pose[first_frame_id].data(),
                                     &pt.second.Idepth_opt);
        }

        // add point constraints
        {
            for (const auto &ob: obs) {
                if (ob.first == first_frame_id) continue;

                ceres::LossFunction *loss_function = nullptr;
                loss_function = new ceres::CauchyLoss(1.);
                auto cost_function = new ProjectionPointInverseDepthFactor(obs.at(first_frame_id).normalpoint,
                                                                           ob.second.normalpoint);
                problem.AddResidualBlock(cost_function, loss_function, para_Pose[first_frame_id].data(),
                                         para_Pose[ob.first].data(), &pt.second.Idepth_opt);

            }
        }
    }

    {
        // plane prior constrain
        for (int i = 0; i < planes_opt.size(); i++) {
            ceres::LossFunction *loss = nullptr;
            loss = new ceres::CauchyLoss(1.);

            ceres::CostFunction* cost_prior_plane = nullptr;
            if(constraint_type == 0) {
                cost_prior_plane = new PointOnPlaneProjectionPriorFactorWall(planes_prior[i]);
            }else
            {
                cost_prior_plane = new PointOnPlaneProjectionPriorFactorGround(planes_prior[i]);
            }

            problem.AddResidualBlock(cost_prior_plane, loss, planes_opt[i].data());
        }
    }

    ceres::Solver::Options options;
//    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;  // LEVENBERG_MARQUARDT  DOGLEG
    // options.linear_solver_type = ceres::SPARSE_SCHUR; // SPARSE_NORMAL_CHOLESKY  or SPARSE_SCHUR
    options.max_num_iterations = optimization_num;
    options.minimizer_progress_to_stdout = false;

    TicToc solver_time;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // std::cout << "loose solver time: " <<  solver_time.toc() << std::endl;
    // info.opti_step = summary.iterations.size();
    std::cout << summary.BriefReport() << std::endl;

    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.block(0, 3, 3, 1) = tic;
    Tic.block(0, 0, 3, 3) = Ric;

    Twis_opt.resize(Twcs.size());
    Twcs_opt.resize(Twcs.size());

    for (int i = 0; i < Twcs.size(); ++i) {
        Eigen::Quaterniond qua(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]);
        Eigen::Matrix3d rot = qua.toRotationMatrix();
        Eigen::Vector3d tran(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);

        Eigen::Matrix4d Twi = Eigen::Matrix4d::Identity();
        Twi.block(0, 3, 3, 1) = tran;
        Twi.block(0, 0, 3, 3) = rot;
        Twcs_opt[i] = Twi * Tic;
        Twis_opt[i] = Twi;
    }


}

void BundleAdjustmentBase::tightly_coupled_method() {

    // weight of two constraints
    ProjectionPointInverseDepthFactor::sqrt_info = 460.0 / 1.5 * Eigen::Matrix2d::Identity();
    PointOnPlaneProjectionFactor::sqrt_info = 460.0 / 0.5 * Eigen::Matrix2d::Identity();
    PointOnPlaneProjectionPriorFactorWall::sqrt_info = 1e5 * Eigen::Matrix3d::Identity();
    PointOnPlaneProjectionPriorFactorGround::sqrt_info = 1e5 * Eigen::Matrix3d::Identity();

    ProjectionPointInverseDepthFactor::qic = Eigen::Quaterniond(Ric);
    ProjectionPointInverseDepthFactor::tic = tic;

    PointOnPlaneProjectionFactor::qic = Eigen::Quaterniond(Ric);
    PointOnPlaneProjectionFactor::tic = tic;

    ceres::Problem problem;

    planes_opt.clear();
    planes_opt = planes;

    // plane parameters
    for (int i = 0; i < planes_opt.size(); i++) {
        ceres::LocalParameterization *local_parameterization;
        local_parameterization = new mesh::WallLocalParameterization();
        problem.AddParameterBlock(planes_opt[i].data(), 4, local_parameterization);

    }

    // pose parameters
    std::vector<Eigen::Matrix<double, 7, 1>> para_Pose(Twis.size());
    for (int i = 0; i < Twis.size(); ++i) {

        Eigen::Matrix4d Twi = Twis[i];
        Eigen::Matrix3d Rot = Twi.block(0, 0, 3, 3);
        Eigen::Vector3d Tran = Twi.block(0, 3, 3, 1);
        Eigen::Quaterniond qua(Rot);

        para_Pose[i][0] = Tran(0);
        para_Pose[i][1] = Tran(1);
        para_Pose[i][2] = Tran(2);
        para_Pose[i][3] = qua.x();
        para_Pose[i][4] = qua.y();
        para_Pose[i][5] = qua.z();
        para_Pose[i][6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i].data(), 7, local_parameterization);

        // fix 10 frames
        if (i < 10)
            problem.SetParameterBlockConstant(para_Pose[i].data());
    }

    for (auto &pt: points_obs) {

        // 根据特征点的点号找到落在的平面位置
        auto plane_id = pt.second.plane_id;
        // 第一帧对应的帧号
        TimeFrameId first_frame_id = pt.second.kf_id;

        auto &obs = pt.second.obs;

        // add plane constraints
        {
            for (const auto &ob: obs) {
                if (ob.first == first_frame_id) continue;
                ceres::LossFunction *loss = nullptr;
                loss = new ceres::CauchyLoss(1.);
                auto cost_plane = new PointOnPlaneProjectionFactor(obs.at(first_frame_id).normalpoint,
                                                                   ob.second.normalpoint);
                problem.AddResidualBlock(cost_plane, loss, planes_opt[plane_id].data(),
                                         para_Pose[first_frame_id].data(),
                                         para_Pose[ob.first].data());
            }
        }

        {
            // plane prior constrain
            for (int i = 0; i < planes_opt.size(); i++) {
                ceres::LossFunction *loss = nullptr;
                loss = new ceres::CauchyLoss(1.);

                ceres::CostFunction* cost_prior_plane = nullptr;
                if(constraint_type == 0) {
                    cost_prior_plane = new PointOnPlaneProjectionPriorFactorWall(planes_prior[i]);
                }else
                {
                    cost_prior_plane = new PointOnPlaneProjectionPriorFactorGround(planes_prior[i]);
                }
                problem.AddResidualBlock(cost_prior_plane, loss, planes_opt[i].data());
            }
        }

    }


    ceres::Solver::Options options;
//    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;  // LEVENBERG_MARQUARDT  DOGLEG
    // options.linear_solver_type = ceres::SPARSE_SCHUR; // SPARSE_NORMAL_CHOLESKY  or SPARSE_SCHUR
    options.max_num_iterations = optimization_num;
    options.minimizer_progress_to_stdout = false;

    TicToc solver_time;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // std::cout << "tight solver time: " <<  solver_time.toc() << std::endl;
    // info.opti_step = summary.iterations.size();
    std::cout << summary.BriefReport() << std::endl;

    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.block(0, 3, 3, 1) = tic;
    Tic.block(0, 0, 3, 3) = Ric;

    Twis_opt.resize(Twcs.size());
    Twcs_opt.resize(Twcs.size());

    for (int i = 0; i < Twcs.size(); ++i) {
        Eigen::Quaterniond qua(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]);
        Eigen::Matrix3d rot = qua.toRotationMatrix();
        Eigen::Vector3d tran(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);

        Eigen::Matrix4d Twi = Eigen::Matrix4d::Identity();
        Twi.block(0, 3, 3, 1) = tran;
        Twi.block(0, 0, 3, 3) = rot;
        Twcs_opt[i] = Twi * Tic;
        Twis_opt[i] = Twi;
    }

    TicToc t_update_depth;
    Eigen::Vector3d optical_center{0., 0., 0.};
    for (auto &pt: points_obs) {
        PlaneID plane_id = pt.second.plane_id;
        TimeFrameId host_id = pt.second.kf_id;
        // nc = Rcw * nw
        Eigen::Vector3d norm_c = Twcs[host_id].block<3, 3>(0, 0).transpose() * planes_opt[plane_id].head<3>();
        // dc = dw + twc * nw
        double distance_c =
                planes_opt[plane_id](3) + Twcs[host_id].block<3, 1>(0, 3).dot(planes_opt[plane_id].head<3>());

        pt.second.Idepth = updateCoplanarPoint(optical_center, pt.second.obs.at(host_id).normalpoint,
                                               Eigen::Vector4d(norm_c.x(), norm_c.y(), norm_c.z(), distance_c));
    }
}


double BundleAdjustmentBase::updateCoplanarPoint(const Eigen::Vector3d &lineStart, const Eigen::Vector3d &lineEnd,
                                                 const Eigen::Vector4d &plane_para_host) {

    Eigen::Vector3d dir = lineEnd - lineStart;
    double b = lineStart.dot(plane_para_host.head<3>()) + (plane_para_host.w());
    double a = dir.dot(plane_para_host.head<3>());
    double t = -b / a;
    Eigen::Vector3d res = lineStart + t * dir;
    return 1 / res.z();
}

