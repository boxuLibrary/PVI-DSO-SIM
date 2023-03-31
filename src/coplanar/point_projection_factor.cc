
#include "point_projection_factor.hpp"
#include "plane_local_parameterization.hpp"
#include "tic_toc.h"


Eigen::Matrix2d ProjectionPointInverseDepthFactor::sqrt_info;
double ProjectionPointInverseDepthFactor::sum_t;
Eigen::Quaterniond ProjectionPointInverseDepthFactor::qic;
Eigen::Vector3d ProjectionPointInverseDepthFactor::tic;


ProjectionPointInverseDepthFactor::ProjectionPointInverseDepthFactor(const Eigen::Vector3d &ob0,
                                                                     const Eigen::Vector3d &obi) : point_obi_(obi),
                                                                                                   point_ob0_(ob0) {}

bool ProjectionPointInverseDepthFactor::Evaluate(double const *const *parameters, double *residuals,
                                                 double **jacobians) const {

    Eigen::Vector3d twi0(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond qwi0(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d twi1(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qwi1(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    double inverse_depth = parameters[2][0];

    Eigen::Vector3d pts_camera_0 = point_ob0_ / inverse_depth;
    Eigen::Vector3d pts_imu_0 = qic * pts_camera_0 + tic;
    Eigen::Vector3d pts_w = qwi0 * pts_imu_0 + twi0;
    Eigen::Vector3d pts_imu_1 = qwi1.inverse() * (pts_w - twi1);
    Eigen::Vector3d pts_camera_1 = qic.inverse() * (pts_imu_1 - tic);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual(0) = pts_camera_1(0) / pts_camera_1(2) - point_obi_(0);
    residual(1) = pts_camera_1(1) / pts_camera_1(2) - point_obi_(1);

    residual = sqrt_info * residual;


    if (jacobians) {

        Eigen::Matrix<double, 2, 3> reduce(2, 3);
        reduce << 1. / pts_camera_1(2), 0, -pts_camera_1(0) / (pts_camera_1(2) * pts_camera_1(2)),
                0, 1. / pts_camera_1(2), -pts_camera_1(1) / (pts_camera_1(2) * pts_camera_1(2));
        reduce = sqrt_info * reduce;

        Eigen::Matrix3d Rwc0 = qwi0.toRotationMatrix() * qic.toRotationMatrix();
        Eigen::Matrix3d Rwc1 = qwi1.toRotationMatrix() * qic.toRotationMatrix();

        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose0(jacobians[0]);
            Eigen::Matrix3d Inde = Eigen::Matrix3d::Identity();
            Eigen::Matrix<double, 3, 6> jaco_pose0;
            jaco_pose0.leftCols<3>() = Rwc1.transpose();
            jaco_pose0.rightCols<3>() = Rwc1.transpose() *
                                        qwi0.toRotationMatrix() * -Utility::skewSymmetric(pts_imu_0);

            jacobian_pose0.leftCols<6>() = reduce * jaco_pose0;
            jacobian_pose0.rightCols<1>().setZero();

        }

        if (jacobians[1]) {
            TicToc t_posej;
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose1(jacobians[1]);
            Eigen::Matrix3d Inde = Eigen::Matrix3d::Identity();
            Eigen::Matrix<double, 3, 6> jaco_pose1;

            jaco_pose1.leftCols<3>() = -Rwc1.transpose();
            jaco_pose1.rightCols<3>() = qic.toRotationMatrix().transpose() * Utility::skewSymmetric(pts_imu_1);

            jacobian_pose1.leftCols<6>() = reduce * jaco_pose1;
            jacobian_pose1.rightCols<1>().setZero();
        }

        if (jacobians[2]) {
            Eigen::Map<Eigen::Vector2d> jacobian_point(jacobians[2]);
            jacobian_point =
                    reduce * Rwc1.transpose() * Rwc0 * point_ob0_ * -1.0 / (inverse_depth * inverse_depth);

        }



        // Jacobian Check
        if (0 && jacobians[0] && jacobians[1] && jacobians[2]) {
            const double eps = 1e-6;
            Eigen::Matrix<double, 2, 13> num_jacobian;
            for (int i = 0; i < 13; ++i) {

                Eigen::Vector3d twi0_check(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond qwi0_check(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
                Eigen::Vector3d twi1_check(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Quaterniond qwi1_check(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
                double inverse_depth_check = parameters[2][0];

                {
                    int b = i % 3, a = i / 3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                    if (a == 0)
                        twi0_check = twi0_check + delta;
                    else if (a == 1)
                        qwi0_check = qwi0_check * Utility::deltaQ(delta);
                    else if (a == 2)
                        twi1_check = twi1_check + delta;
                    else if (a == 3)
                        qwi1_check = qwi1_check * Utility::deltaQ(delta);
                    else
                        inverse_depth_check += eps;
                }

                Eigen::Vector3d twc0_check = twi0_check + qwi0_check * tic;
                Eigen::Quaterniond qwc0_check = qwi0_check * qic;
                Eigen::Vector3d twc1_check = twi1_check + qwi1_check * tic;
                Eigen::Quaterniond qwc1_check = qwi1_check * qic;

                Eigen::Vector3d pts_camera_0_check = point_ob0_ / inverse_depth_check;
                Eigen::Vector3d pts_w_check = qwc0_check * pts_camera_0_check + twc0_check;
                Eigen::Vector3d pts_camera_i_check = qwc1_check.inverse() * (pts_w_check - twc1_check);

                Eigen::Vector2d residual_tmp;
                residual_tmp(0) = pts_camera_i_check(0) / pts_camera_i_check(2) - point_obi_(0);
                residual_tmp(1) = pts_camera_i_check(1) / pts_camera_i_check(2) - point_obi_(1);

                residual_tmp = sqrt_info * residual_tmp;
                num_jacobian.col(i) = (residual_tmp - residual) / eps;
            }

            // check jacobian
            std::cout << "ana = " << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[0]) << std::endl
                      << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[1]) << std::endl
                      << std::endl;
            std::cout << Eigen::Map<Eigen::Vector2d>(jacobians[2]) << std::endl
                      << std::endl;
            std::cout << "num_jacobian:\n" << num_jacobian.leftCols(6) << "\n" << num_jacobian.rightCols(7) << "\n"
                      << std::endl;

            std::cout << " ------------------- " << std::endl;
        }
    }

    return true;
}


double PointOnPlaneFactor::sqrt_info;// = 1.0;
Eigen::Quaterniond PointOnPlaneFactor::qic;
Eigen::Vector3d PointOnPlaneFactor::tic;

PointOnPlaneFactor::PointOnPlaneFactor(const Eigen::Vector3d &_pts_i) : pts_i(_pts_i) {

}

/*
  parameters[0]:  Twi
  parameters[1]:  Twj
  parameters[2]:  Tbc
  parameters[3]:  line_orth
*/
bool PointOnPlaneFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {

    double distance = parameters[0][3];
    Eigen::Vector3d nom_n_now(parameters[0][0], parameters[0][1], parameters[0][2]);

    Eigen::Vector3d Pi(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qi(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    double inv_dep_i = parameters[2][0];

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;

    residuals[0] = sqrt_info * (nom_n_now.dot(pts_w) + distance);

    if (jacobians) {
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 1, 4, Eigen::RowMajor>> jacobian_plane(jacobians[0]);
            jacobian_plane(0, 0) = pts_w.x();
            jacobian_plane(0, 1) = pts_w.y();
            jacobian_plane(0, 2) = pts_w.z();
            jacobian_plane(0, 3) = 1;
            jacobian_plane = sqrt_info * jacobian_plane;
        }


        Eigen::Vector3d jacobian_word(nom_n_now(0), nom_n_now(1), nom_n_now(2));

        jacobian_word = sqrt_info * jacobian_word;

        if (jacobians[1]) {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[1]);

            jacobian_pose.leftCols<3>() = jacobian_word.transpose();

            Eigen::Matrix3d Qi_R = Qi.toRotationMatrix();
            Eigen::Matrix<double, 1, 3> jacob_Qi =
                    -jacobian_word.transpose() * Qi_R * Utility::skewSymmetric(pts_imu_i);
            jacobian_pose(0, 3) = jacob_Qi(0);
            jacobian_pose(0, 4) = jacob_Qi(1);
            jacobian_pose(0, 5) = jacob_Qi(2);
            jacobian_pose(0, 6) = 0;
        }


        Eigen::Matrix3d Qi_R = Qi.toRotationMatrix();
        Eigen::Matrix<double, 1, 3> jacobian_imu = jacobian_word.transpose() * Qi_R;

        Eigen::Matrix3d qic_R = qic.toRotationMatrix();
        Eigen::Matrix<double, 1, 3> jacobian_camera = jacobian_imu * qic_R;

        if (jacobians[2]) {
            double jac = jacobian_camera * pts_i;
            jacobians[2][0] = -jac / (inv_dep_i * inv_dep_i);
        }

        // Jacobian Check
        if (0 && jacobians[0] && jacobians[1] && jacobians[2]) {
            const double eps = 1e-9;
            Eigen::Matrix<double, 1, 17> num_jacobian;
            for (int i = 0; i < 17; ++i) {

                Eigen::Vector3d plane_check(parameters[0][0], parameters[0][1], parameters[0][2]);

                double distance_check = parameters[0][3];

                Eigen::Vector3d Pi_check(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Quaterniond Qi_check(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                double inv_dep_i_check = parameters[2][0];

                Eigen::Vector3d tic_check = tic;
                Eigen::Quaterniond qic_check = qic;


                if (i > 3) {
                    int a = (i - 4) / 3, b = (i - 4) % 3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                    if (a == 0)
                        Pi_check += delta;
                    else if (a == 1)
                        Qi_check = Qi_check * Utility::deltaQ(delta);
                    else if (a == 2)
                        tic_check += delta;
                    else if (a == 3)
                        qic_check = qic_check * Utility::deltaQ(delta);
                    else if (a == 4)
                        inv_dep_i_check += delta.x();
                } else {
                    int a = i / 3, b = i % 3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

                    if (a == 0)
                        plane_check += delta;
                    else if (a == 1)
                        distance_check += delta.x();

                }


                Eigen::Vector3d nom_n_now_check = plane_check;
                Eigen::Vector3d pts_camera_i_check = pts_i / inv_dep_i_check;
                Eigen::Vector3d pts_imu_i_check = qic_check * pts_camera_i_check + tic_check;
                Eigen::Vector3d pts_w_check = Qi_check * pts_imu_i_check + Pi_check;

                // 误差
                double tmp_residual = sqrt_info * (nom_n_now_check.dot(pts_w_check) + distance_check);

                num_jacobian[i] = (tmp_residual - residuals[0]) / eps;

                std::cout << " point : error " << residuals[0] << ":" << tmp_residual << std::endl;
            }

            std::cout << "-----analysis------:" << std::endl;
            std::cout << "plane: ";
            for (int i = 0; i < 4; ++i)
                std::cout << jacobians[0][i] << "  ";
            std::cout << std::endl;
            std::cout << "pose: ";
            for (int i = 0; i < 6; ++i)
                std::cout << jacobians[1][i] << "  ";
            std::cout << std::endl;
            std::cout << "inverse depth: " << jacobians[2][0] << std::endl;

            std::cout << "-------num--------:" << std::endl;
            std::cout << "plane: ";
            for (int i = 0; i < 4; ++i)
                std::cout << num_jacobian[i] << "  ";
            std::cout << std::endl;
            std::cout << "pose: ";
            for (int i = 4; i < 10; ++i)
                std::cout << num_jacobian[i] << "  ";
            std::cout << std::endl;
            std::cout << "expose: ";
            for (int i = 10; i < 16; ++i)
                std::cout << num_jacobian[i] << "  ";
            std::cout << std::endl;
            std::cout << "inverse depth: ";
            std::cout << num_jacobian[16] << "  " << std::endl;
            std::cout << std::endl;
            std::cout << " ------------------- " << std::endl;

        }

    }

    return true;
}


Eigen::Matrix2d PointOnPlaneProjectionFactor::sqrt_info;// = 1.0;
Eigen::Quaterniond PointOnPlaneProjectionFactor::qic;// = 1.0;
Eigen::Vector3d PointOnPlaneProjectionFactor::tic;// = 1.0;

PointOnPlaneProjectionFactor::PointOnPlaneProjectionFactor(
        const Eigen::Vector3d &_ob0,
        const Eigen::Vector3d &_obi) : ob0(_ob0), obi(_obi) {}


/*
  parameters[0]:  Twi
  parameters[1]:  Twj
  parameters[2]:  Tbc
  parameters[3]:  line_orth
*/
bool
PointOnPlaneProjectionFactor::Evaluate(double const *const *parameters, double *residuals,
                                       double **jacobians) const {

    // TicToc t_solve1;
    // plane
    Eigen::Vector3d norm_w(parameters[0][0], parameters[0][1], parameters[0][2]);
    double distance_w = parameters[0][3];
    // pose0
    Eigen::Vector3d Pi0(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qi0(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    // pose1
    Eigen::Vector3d Pi1(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond Qi1(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    const Eigen::Quaterniond q_w_c_h = Qi0 * qic;
    const Eigen::Quaterniond q_w_c_t = Qi1 * qic;
    const Eigen::Vector3d &t_w_c_h = Pi0 + Qi0 * tic;
    const Eigen::Vector3d &t_w_c_t = Pi1 + Qi1 * tic;
    const Eigen::Quaterniond q_t_h = q_w_c_t.conjugate() * q_w_c_h;

    // nc = Rcw * nw
    Eigen::Vector3d norm_c = q_w_c_h.conjugate() * norm_w;
    // dc = dw + twc * nw
    double distance_c = distance_w + t_w_c_h.dot(norm_w);
    double depth_cam = -distance_c / (norm_c.dot(ob0));

    Eigen::Vector3d pts_camera_0 = ob0 * depth_cam;
    Eigen::Vector3d pts_w = q_w_c_h * pts_camera_0 + t_w_c_h;
    Eigen::Vector3d pts_camera_i = q_w_c_t.conjugate() * (pts_w - t_w_c_t);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual(0) = pts_camera_i(0) / pts_camera_i(2) - obi(0);
    residual(1) = pts_camera_i(1) / pts_camera_i(2) - obi(1);
    residual = sqrt_info * residual;

    if (jacobians) {

        //  TicToc t_solve2;
        // d_uv_xyz
        Eigen::Matrix<double, 2, 3> jaco_xyz;
        jaco_xyz << 1 / pts_camera_i.z(), 0, -pts_camera_i.x() / (pts_camera_i.z() * pts_camera_i.z()),
                0, 1 / pts_camera_i.z(), -pts_camera_i.y() / (pts_camera_i.z() * pts_camera_i.z());

        jaco_xyz = sqrt_info * jaco_xyz;

        Eigen::Vector3d nw = norm_w;
        double dw = distance_w;
        double n_w_R_w_i_pts = norm_c.dot(ob0);


        Eigen::Vector3d ob_t = q_t_h * ob0;
        Eigen::Vector3d ob_w = q_w_c_h * ob0;
        Eigen::Vector3d ob_imu = qic * ob0;

        Eigen::Matrix3d R_c_w_t = q_w_c_t.conjugate().toRotationMatrix();


        if (jacobians[0]) {
            // TicToc t_plane;

            Eigen::Map<Eigen::Matrix<double, 2, 4, Eigen::RowMajor>> jacobian_plane(jacobians[0]);
            // plane jacobian
            // d_xyz_nd
            Eigen::Matrix3d S1 = ob_t * t_w_c_h.transpose();
            Eigen::Matrix<double, 1, 3> S2 = ob_w.transpose();

            Eigen::Matrix3d jaco_global_normal =
                    -(S1 * n_w_R_w_i_pts - ob_t * distance_c * S2) /
                    (n_w_R_w_i_pts * n_w_R_w_i_pts);

            Eigen::Matrix<double, 3, 4> jaco_plane_w;
            jaco_plane_w.block<3, 3>(0, 0) = jaco_global_normal;
            jaco_plane_w.col(3) = -ob_t / n_w_R_w_i_pts;
            jacobian_plane = jaco_xyz * jaco_plane_w;

        }

        if (jacobians[1]) {

            // TicToc t_posei;

            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posei(jacobians[1]);

            // 1. T_w_c_h translation;
            Eigen::Matrix3d jaco_host_trans =
                    -ob_t * nw.transpose() / n_w_R_w_i_pts + R_c_w_t;

            // 2. T_w_c_h rotation
            Eigen::Matrix3d A1 = -R_c_w_t * Qi0.toRotationMatrix() * Utility::skewSymmetric(ob_imu * distance_c);

            Eigen::Matrix<double, 1, 3> A2 = ob_imu.transpose() * Utility::skewSymmetric(Qi0.conjugate() * norm_w);

            Eigen::Matrix<double, 3, 1> A3 = distance_c * ob_t;

            Eigen::Matrix3d jaco_host_rot = -(A1 * n_w_R_w_i_pts - A3 * A2) / (n_w_R_w_i_pts * n_w_R_w_i_pts);

            jacobian_posei.block<2, 3>(0, 0) = jaco_xyz * jaco_host_trans;
            jacobian_posei.block<2, 3>(0, 3) = jaco_xyz * jaco_host_rot;

            jacobian_posei.rightCols<1>().setZero();

        }

        if (jacobians[2]) {

            // TicToc t_posej;

            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posej(jacobians[2]);
            // 3. T_w_c_t translation
            Eigen::Matrix3d jaco_target_trans = -R_c_w_t;
            // 4. T_w_c_t rotation
            Eigen::Vector3d jaRj_tmp;
            jaRj_tmp = -ob_w * distance_c / n_w_R_w_i_pts;
            jaRj_tmp = jaRj_tmp + t_w_c_h - t_w_c_t;

            Eigen::Vector3d M = Qi1.conjugate().toRotationMatrix() * jaRj_tmp;

            Eigen::Matrix3d jaco_target_rot = qic.conjugate().toRotationMatrix() * Utility::skewSymmetric(M);

            jacobian_posej.block<2, 3>(0, 0) = jaco_xyz * jaco_target_trans;
            jacobian_posej.block<2, 3>(0, 3) = jaco_xyz * jaco_target_rot;
            jacobian_posej.rightCols<1>().setZero();


        }

        // Jacobian Check
        if (0 && jacobians[0] && jacobians[1] && jacobians[2]) {

            const double eps = 1e-6;
            Eigen::Matrix<double, 2, 22> num_jacobian;

            TicToc t_solver3;
            for (int i = 0; i < 22; ++i) {

                double distance_w_ck = parameters[0][3];
                Eigen::Vector3d normal_w_ck(parameters[0][0], parameters[0][1], parameters[0][2]);

                Eigen::Vector3d Pi0_ck(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Quaterniond Qi0_ck(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                Eigen::Vector3d Pi1_ck(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Quaterniond Qi1_ck(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

                Eigen::Vector3d tic_ck = tic;
                Eigen::Quaterniond qic_ck = qic;

                if (i > 3) {
                    int b = (i - 4) % 3, a = (i - 4) / 3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                    if (a == 0) {
                        Pi0_ck = Pi0_ck + delta;
                    } else if (a == 1)
                        Qi0_ck = Qi0_ck * Utility::deltaQ(delta);
                    else if (a == 2)
                        Pi1_ck = Pi1_ck + delta;
                    else if (a == 3)
                        Qi1_ck = Qi1_ck * Utility::deltaQ(delta);
                    else if (a == 4)
                        tic_ck = tic_ck + delta;
                    else if (a == 5)
                        qic_ck = qic_ck * Utility::deltaQ(delta);
                } else {

                    int b = i % 3, a = i / 3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                    if (a == 0) {
                        normal_w_ck = normal_w_ck + delta;
                    } else if (a == 1)
                        distance_w_ck += delta.x();

                }

                Sophus::SE3d T_i_c_ck(qic_ck, tic_ck);
                Sophus::SE3d T_w_i_h_ck(Qi0_ck, Pi0_ck);
                Sophus::SE3d T_w_i_t_ck(Qi1_ck, Pi1_ck);

                Sophus::SE3d T_w_c_h_ck = T_w_i_h_ck * T_i_c_ck;
                Sophus::SE3d T_w_c_t_ck = T_w_i_t_ck * T_i_c_ck;
                Sophus::SE3d T_c_t_c_h_ck = T_w_c_t_ck.inverse() * T_w_c_h_ck;

                Eigen::Vector3d norm_c_ck = T_w_c_h_ck.inverse().rotationMatrix() * normal_w_ck;

                double distance_c_ck = distance_w_ck + T_w_c_h_ck.translation().dot(normal_w_ck);

                double depth_cam_ck = -distance_c_ck / (norm_c_ck.dot(ob0));

                Eigen::Vector3d pts_camera_0_ck = ob0 * depth_cam_ck;
                Eigen::Vector3d pts_w_ck = T_w_c_h_ck.rotationMatrix() * pts_camera_0_ck + T_w_c_h_ck.translation();
                Eigen::Vector3d pts_camera_i_ck =
                        T_w_c_t_ck.rotationMatrix().inverse() * (pts_w_ck - T_w_c_t_ck.translation());

                Eigen::Vector2d residual_tmp;
                residual_tmp(0) = pts_camera_i_ck(0) / pts_camera_i_ck(2) - obi(0);
                residual_tmp(1) = pts_camera_i_ck(1) / pts_camera_i_ck(2) - obi(1);

                residual_tmp = sqrt_info * residual_tmp;
                num_jacobian.col(i) = (residual_tmp - residual) / eps;


//                const Eigen::Quaterniond q_w_c_h_ck = Qi0_ck * qic_ck;
//                const Eigen::Quaterniond q_w_c_t_ck = Qi1_ck * qic_ck;
//                const Eigen::Vector3d &t_w_c_h_ck = Pi0_ck + Qi0_ck * tic_ck;
//                const Eigen::Vector3d &t_w_c_t_ck = Pi1_ck + Qi1_ck * tic_ck;
//                const Eigen::Quaterniond q_t_h_ck = q_w_c_t_ck.conjugate() * q_w_c_h_ck;
//
//                // nc = Rcw * nw
//                Eigen::Vector3d norm_c_ck = q_w_c_h_ck.conjugate() * normal_w_ck;
//                // dc = dw + twc * nw
//                double distance_c_ck = distance_w_ck + t_w_c_h_ck.dot(normal_w_ck);
//                double depth_cam_ck = -distance_c_ck / (norm_c_ck.dot(ob0));
//
//                Eigen::Vector3d pts_camera_0_ck = ob0 * depth_cam_ck;
//                Eigen::Vector3d pts_w_ck = q_w_c_h_ck * pts_camera_0_ck + t_w_c_h_ck;
//                Eigen::Vector3d pts_camera_i_ck = q_w_c_t_ck.conjugate() * (pts_w_ck - t_w_c_t_ck);
//
//
//                Eigen::Vector2d residual_tmp;
//                residual_tmp(0) = pts_camera_i_ck(0) / pts_camera_i_ck(2) - obi(0);
//                residual_tmp(1) = pts_camera_i_ck(1) / pts_camera_i_ck(2) - obi(1);
//
//                residual_tmp = sqrt_info * residual_tmp;
//                num_jacobian.col(i) = (residual_tmp - residual) / eps;
            }

            std::cout << "t_solve3: " << t_solver3.toc() << std::endl;
            // check jacobian
            std::cout << "ana = " << std::endl;
            std::cout << "plane" << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 4, Eigen::RowMajor>>(jacobians[0]) << std::endl;
            std::cout << "pose0" << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[1]) << std::endl;
            std::cout << "pose1" << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[2]) << std::endl;


            std::cout << "num_jacobian:\n";
            std::cout << "plane" << std::endl;
            std::cout << num_jacobian.block<2, 4>(0, 0) << "\n";
            std::cout << "pose0" << std::endl;
            std::cout << num_jacobian.block<2, 6>(0, 4) << "\n";
            std::cout << "pose1" << std::endl;
            std::cout << num_jacobian.block<2, 6>(0, 10) << "\n";
            std::cout << "expose: " << std::endl;
            std::cout << num_jacobian.block<2, 6>(0, 16) << std::endl;
            std::cout << " ------------------- " << std::endl;
        }
    }

    return true;
}


Eigen::Matrix3d PointOnPlaneProjectionPriorFactorWall::sqrt_info;// = 1.0;


PointOnPlaneProjectionPriorFactorWall::PointOnPlaneProjectionPriorFactorWall(
        const Eigen::Vector4d &plane_prior_) : plane_prior(plane_prior_) {}


/*
  parameters[0]:  Twi
  parameters[1]:  Twj
  parameters[2]:  Tbc
  parameters[3]:  line_orth
*/
bool
PointOnPlaneProjectionPriorFactorWall::Evaluate(double const *const *parameters, double *residuals,
                                            double **jacobians) const {

    // TicToc t_solve1;
    // plane

    Eigen::Vector3d norm_w(parameters[0][0], parameters[0][1], parameters[0][2]);
    double distance_w = parameters[0][3];

    double theta = std::atan2(norm_w.y(), norm_w.x());
    double phi = std::asin(norm_w.z());

    Eigen::Vector3d norm_obs = plane_prior.head(3);
    double distance_obs = plane_prior(3);

    double theta_obs = std::atan2(norm_obs.y(), norm_obs.x());
    double phi_obs = std::asin(norm_obs.z());

    Eigen::Map<Eigen::Vector3d> residual(residuals);
    residual(0) = theta - theta_obs;
    residual(1) = phi - phi_obs;
    residual(2) = distance_w - distance_obs;

    residual = sqrt_info * residual;


    if (jacobians) {

        if (jacobians[0]) {
            // TicToc t_plane;
            Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> jacobian_plane(jacobians[0]);

            double nx = norm_w.x();
            double ny = norm_w.y();
            double nz = norm_w.z();

            double base1 = std::pow(nx, 2) / (std::pow(nx, 2) + std::pow(ny, 2));
            double base2 = 1. / std::sqrt(1 - std::pow(nz, 2));

            Eigen::Matrix<double, 4, 1> jaco0, jaco1, jaco2;

            jaco0 << base1 * (-ny / std::pow(nx, 2)), base1 * (1. / nx), 0., 0.;
            jaco1 << 0., 0., base2, 0.;
            jaco2 << 0., 0., 0., 1.;

            jacobian_plane.row(0) = jaco0;
            jacobian_plane.row(1) = jaco1;
            jacobian_plane.row(2) = jaco2;
            jacobian_plane = sqrt_info * jacobian_plane;

        }


        // Jacobian Check
        if (0 && jacobians[0]) {

            const double eps = 1e-6;
            Eigen::Matrix<double, 3, 4> num_jacobian;

            for (int i = 0; i < 4; ++i) {

                double distance_w_ck = parameters[0][3];
                Eigen::Vector3d normal_w_ck(parameters[0][0], parameters[0][1], parameters[0][2]);


                int a = i / 3, b = i % 3;
                Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                if (a == 0) {
                    normal_w_ck += delta;
                } else if (a == 1) {
                    distance_w_ck += delta.x();
                }


                double theta_ck = std::atan2(normal_w_ck.y(), normal_w_ck.x());
                double phi_ck = std::asin(normal_w_ck.z());


                Eigen::Vector3d residual_tmp;
                residual_tmp(0) = theta_ck - theta_obs;
                residual_tmp(1) = phi_ck - phi_obs;
                residual_tmp(2) = distance_w_ck - distance_obs;

                residual_tmp = sqrt_info * residual_tmp;
                num_jacobian.col(i) = (residual_tmp - residual) / eps;

            }

            // check jacobian

            std::cout << "ana = " << std::endl;
            std::cout << "plane prior" << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>>(jacobians[0]) << std::endl;

            std::cout << "num_jacobian:\n";
            std::cout << "plane prior" << std::endl;
            std::cout << num_jacobian.block<3, 4>(0, 0) << "\n";


        }
    }

    return true;
}


Eigen::Matrix3d PointOnPlaneProjectionPriorFactorGround::sqrt_info;// = 1.0;


PointOnPlaneProjectionPriorFactorGround::PointOnPlaneProjectionPriorFactorGround(
        const Eigen::Vector4d &plane_prior_) : plane_prior(plane_prior_) {}


/*
  parameters[0]:  Twi
  parameters[1]:  Twj
  parameters[2]:  Tbc
  parameters[3]:  line_orth
*/
bool
PointOnPlaneProjectionPriorFactorGround::Evaluate(double const *const *parameters, double *residuals,
                                            double **jacobians) const {

    // TicToc t_solve1;
    // plane

    // to avoid singular
    Eigen::AngleAxisd rotZ(20. / 180. * M_PI, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd rotX(20. / 180. * M_PI, Eigen::Vector3d::UnitX());


    Eigen::Vector3d norm_w(parameters[0][0], parameters[0][1], parameters[0][2]);
    double distance_w = parameters[0][3];
    norm_w = rotZ * rotX * norm_w;

    double theta = std::atan2(norm_w.y(), norm_w.x());
    double phi = std::asin(norm_w.z());


    Eigen::Vector3d norm_obs = plane_prior.head(3);
    double distance_obs = plane_prior(3);
    norm_obs =  rotZ * rotX * norm_obs;

    double theta_obs = std::atan2(norm_obs.y(), norm_obs.x());
    double phi_obs = std::asin(norm_obs.z());

    Eigen::Map<Eigen::Vector3d> residual(residuals);
    residual(0) = theta - theta_obs;
    residual(1) = phi - phi_obs;
    residual(2) = distance_w - distance_obs;

    residual = sqrt_info * residual;

    if (jacobians) {

        if (jacobians[0]) {
            // TicToc t_plane;
            Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> jacobian_plane(jacobians[0]);

            double nx = norm_w.x();
            double ny = norm_w.y();
            double nz = norm_w.z();

            double base1 = std::pow(nx, 2) / (std::pow(nx, 2) + std::pow(ny, 2));
            double base2 = 1. / std::sqrt(1 - std::pow(nz, 2));

            Eigen::Matrix<double, 4, 1> jaco0, jaco1, jaco2;

            jaco0 << base1 * (-ny / std::pow(nx, 2)), base1 * (1. / nx), 0., 0.;
            jaco1 << 0., 0., base2, 0.;
            jaco2 << 0., 0., 0., 1.;

            jacobian_plane.row(0) = jaco0;
            jacobian_plane.row(1) = jaco1;
            jacobian_plane.row(2) = jaco2;
            jacobian_plane = sqrt_info * jacobian_plane;

        }


        // Jacobian Check
        if (0 && jacobians[0]) {

            const double eps = 1e-6;
            Eigen::Matrix<double, 3, 4> num_jacobian;

            for (int i = 0; i < 4; ++i) {

                double distance_w_ck = parameters[0][3];
                Eigen::Vector3d normal_w_ck(parameters[0][0], parameters[0][1], parameters[0][2]);


                int a = i / 3, b = i % 3;
                Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                if (a == 0) {
                    normal_w_ck += delta;
                } else if (a == 1) {
                    distance_w_ck += delta.x();
                }


                double theta_ck = std::atan2(normal_w_ck.y(), normal_w_ck.x());
                double phi_ck = std::asin(normal_w_ck.z());


                Eigen::Vector3d residual_tmp;
                residual_tmp(0) = theta_ck - theta_obs;
                residual_tmp(1) = phi_ck - phi_obs;
                residual_tmp(2) = distance_w_ck - distance_obs;

                residual_tmp = sqrt_info * residual_tmp;
                num_jacobian.col(i) = (residual_tmp - residual) / eps;

            }

            // check jacobian

            std::cout << "ana = " << std::endl;
            std::cout << "plane prior" << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>>(jacobians[0]) << std::endl;

            std::cout << "num_jacobian:\n";
            std::cout << "plane prior" << std::endl;
            std::cout << num_jacobian.block<3, 4>(0, 0) << "\n";


        }
    }

    return true;
}







