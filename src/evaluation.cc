
#include "evaluation.hpp"

// 这里都是直接用真实值和估计值在对齐轨迹之后，直接比较的
void evaluateTrajectory(std::vector<Eigen::Matrix4d> &Twcs_gt, std::vector<Eigen::Matrix4d> &Twcs_esti) {

    // 位置精度
    Eigen::Matrix<double, 3, Eigen::Dynamic> est_aligned_pose(3, Twcs_esti.size());
    Eigen::Matrix<double, 3, Eigen::Dynamic> gt_aligned_pose(3, Twcs_gt.size());

    for (int i = 0; i < Twcs_gt.size(); i++) {
        est_aligned_pose(0, i) = Twcs_esti[i].block<3, 1>(0, 3)(0);
        est_aligned_pose(1, i) = Twcs_esti[i].block<3, 1>(0, 3)(1);
        est_aligned_pose(2, i) = Twcs_esti[i].block<3, 1>(0, 3)(2);

        gt_aligned_pose(0, i) = Twcs_gt[i].block<3, 1>(0, 3)(0);
        gt_aligned_pose(1, i) = Twcs_gt[i].block<3, 1>(0, 3)(1);
        gt_aligned_pose(2, i) = Twcs_gt[i].block<3, 1>(0, 3)(2);
    }

    Eigen::Matrix4d Tts = Eigen::umeyama(est_aligned_pose, gt_aligned_pose, true);
    Eigen::Matrix3d cR = Tts.block<3, 3>(0, 0);
    Eigen::Vector3d t = Tts.block<3, 1>(0, 3);
    double s = cR.determinant();
    s = pow(s, 1.0 / 3);
    Eigen::Matrix3d R = cR / s;

    // std::cout << "scale is: " << s << std::endl;

    double pose_rmse = 0.;
    for (int i = 0; i < Twcs_gt.size(); i++) {
        Eigen::Vector3d target_pose = R * est_aligned_pose.col(i) + t;

        pose_rmse += (target_pose - gt_aligned_pose.col(i)).dot(target_pose - gt_aligned_pose.col(i));

    }
    pose_rmse /= Twcs_gt.size();
    pose_rmse = std::sqrt(pose_rmse);

    std::cout << "pose rmse: " << pose_rmse << std::endl;

    // 平移误差和旋转误差
    double rot_rmse = 0.;
    for (int i = 0; i < Twcs_gt.size(); i++) {

        Eigen::Matrix3d rij_est = R * Twcs_esti[i].block<3, 3>(0, 0);
        Eigen::Matrix3d rij_gt = Twcs_gt[i].block<3, 3>(0, 0);
        Eigen::Quaterniond qij_est = Eigen::Quaterniond(rij_est);
        Eigen::Quaterniond qij_gt = Eigen::Quaterniond(rij_gt);
        double error =
                std::acos(((qij_gt * qij_est.inverse()).toRotationMatrix().trace() - 1.0) / 2.0) * 180.0 / M_PI;
        rot_rmse += error * error;
    }

    rot_rmse /= Twcs_gt.size();
    rot_rmse = std::sqrt(rot_rmse);

    // std::cout << "rot rmse: " << rot_rmse << std::endl;
}


void evaluatePlane(std::vector<Eigen::Vector4d>& plane_gt, std::vector<Eigen::Vector4d>& plane_est)
{
    double distance_rmse = 0.;
    double normal_rmse = 0.;
    for(int i = 0; i < plane_gt.size(); i++)
    {
        Eigen::Vector4d param_gt = plane_gt[i];
        Eigen::Vector4d param_est = plane_est[i];

        Eigen::Vector3d normal_gt = param_gt.head(3);
        Eigen::Vector3d normal_est = param_est.head(3);

        // std::cout << "normal est: " <<  normal_est.transpose() << std::endl;
        // std::cout << "normal gt: " <<  normal_gt.transpose() << std::endl;

        double normal_error = normal_est.dot(normal_gt) / (normal_est.norm() * normal_gt.norm());

        // std::cout << "normal error: " << normal_error << std::endl;

        normal_error = std::acos(fminl(fmaxl(normal_error,-1.0),1.0)) / M_PI * 180;

        // std::cout << "normal error: " << normal_error << std::endl;

        normal_rmse += normal_error * normal_error;

        double distance_error = param_est(3) - param_gt(3);
        distance_rmse += distance_error * distance_error;
    }

    distance_rmse /= plane_gt.size();
    distance_rmse = std::sqrt(distance_rmse);
    std::cout << "plane distance rmse: " << distance_rmse << std::endl;

    normal_rmse /= plane_gt.size();
    normal_rmse = std::sqrt(normal_rmse);
    std::cout << "plane normal rmse: " << normal_rmse << std::endl;




}