//
// Created by hyj on 17-6-22.
//

#include <fstream>
#include <random>
#include <sys/stat.h>
#include <unordered_map>
#include "imu.h"
#include "viewer.h"
#include "utilities.h"
#include "feature.h"
#include "evaluation.hpp"
#include "tic_toc.h"

#include "linearCoplanarError.h"

using namespace Sophus;

// find 3d planes from a set of 3d points, using sequential ransac
// input: pts,
void find3dPlanes_pts(Points pts, std::map<int, int> &groups,
                      std::vector<Eigen::Vector4d> &planeParams) {
    int maxIterNo = 500;
    int minSolSetSize = 3;
    double pt2planeDistThresh = 1e-3;
    int planeSetSizeThresh = 30;

    // ----- sequential ransac -----
    for (int seq_i = 0; seq_i < 4; ++seq_i) {
        // ---- ransac ----

        if (pts.empty()) break;

        std::vector<int> indexes(pts.size());
        for (int i = 0; i < indexes.size(); ++i) indexes[i] = i;
        std::vector<int> maxInlierSet;
        double max_d;
        Eigen::Vector3d max_n;

        for (int iter = 0; iter < maxIterNo; iter++) {
            std::vector<int> inlierSet;
            random_unique(indexes.begin(), indexes.end(), minSolSetSize); // shuffle

            Eigen::Vector3d base1 = pts[indexes[0]].head(3) - pts[indexes[1]].head(3);
            Eigen::Vector3d base2 = pts[indexes[0]].head(3) - pts[indexes[2]].head(3);

            Eigen::Vector3d n = base1.cross(base2);

            double d = -n.dot(pts[indexes[0]].head(3)); // plane=[n' d];

            if (d > 0) {
                d = -d;
                n = -n;
            }

            for (int i = 0; i < pts.size(); ++i) {

                // condition: dist < pt2planeDistThresh
                double dist = std::abs(n.dot(pts[i].head(3)) + d) / n.norm();
                // double angle = n.dot(pts[i].head(3)) / (n.norm() * pts[i].head(3).norm());
                // condition1: std::abs(acos(angle) / M_PI * 180) <= 50
                if (dist < pt2planeDistThresh) {
                    inlierSet.push_back(i);
                }
            }

            if (inlierSet.size() > maxInlierSet.size()) {
                maxInlierSet = inlierSet;
                max_n = n / n.norm();
                max_d = d / n.norm();
            }

        }

        // std::cout << "maxInlierSet: " << maxInlierSet.size() << std::endl;
        // std::cout << "maxd: " << max_d << " maxn: " << max_n.transpose() << std::endl;

        if (maxInlierSet.size() > planeSetSizeThresh)  // found a new plane
        {
            std::vector<int> plane_associate;// contains gid of coplanar pts
            Eigen::MatrixXd cpPts{4, maxInlierSet.size()};

            for (int i = maxInlierSet.size() - 1; i >= 0; --i) {
                int pt_id = pts[maxInlierSet[i]](3);
                // std::cout << "pt id: " << pt_id << " " << pts[maxInlierSet[i]].head(3).transpose() << std::endl;
                groups[pt_id] = seq_i;
                cpPts.col(i) = Eigen::Vector4d(pts[maxInlierSet[i]].x(), pts[maxInlierSet[i]].y(),
                                               pts[maxInlierSet[i]].z(), 1);
                pts.erase(pts.begin() + maxInlierSet[i]);
            }

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(cpPts.transpose(), Eigen::ComputeThinV);
            Eigen::Vector4d planeParam = svd.matrixV().rightCols<1>();

            planeParam /= planeParam.head(3).norm();

            // ensure the direction of normal vector is consistency with direction of point
            // n * pt + d = 0
            if (planeParam(3) > 0) planeParam = -planeParam;

            planeParams.emplace_back(planeParam);
        }
    }
}

/**
 * @brief Generate 4*n map points from a cube area
 * 
 * @param points 
 * @param d the dist from the center to the side of the cube  
 * @param n number of mappoint 4 * n
 */
void
CreatePoints_wall(Points &points, std::map<int, int> &plane_associate, std::vector<Eigen::Vector4d> &plane_params,
                  double d,
                  int n) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> x_rand(-d, d);
    std::uniform_real_distribution<double> y_rand(-2, 2);
//    std::uniform_real_distribution<double> d_rand(-d/10., d/10.);
    std::uniform_real_distribution<double> d_rand(0., 0.);

    n = std::max(n, 20);
    int id = 1;
    // side 1
    for (size_t i = 0; i < n; i++) {
        Eigen::Vector4d pt_tmp(x_rand(generator), d + d_rand(generator), y_rand(generator), id++);
        points.push_back(pt_tmp);
    }

    // side 2
    for (size_t i = 0; i < n; i++) {
        Eigen::Vector4d pt_tmp(-d + d_rand(generator), x_rand(generator), y_rand(generator), id++);
        points.push_back(pt_tmp);
    }

    // side 3
    for (size_t i = 0; i < n; i++) {
        Eigen::Vector4d pt_tmp(x_rand(generator),  -d + d_rand(generator), y_rand(generator), id++);
        points.push_back(pt_tmp);

    }

    // side 4
    for (size_t i = 0; i < n; i++) {
        Eigen::Vector4d pt_tmp(d + d_rand(generator), x_rand(generator), y_rand(generator), id++);
        points.push_back(pt_tmp);
    }

    Eigen::AngleAxisd rot(30. / 180. * M_PI, Eigen::Vector3d::UnitZ());

    for (auto &pt: points) {
        pt.head<3>() =  rot.toRotationMatrix() * pt.head<3>();
    }


    find3dPlanes_pts(points, plane_associate, plane_params);

    // save points
    save_points("all_points.txt", points);

}


void
CreatePoints_ground(Points &points, std::map<int, int> &plane_associate, std::vector<Eigen::Vector4d> &plane_params,
                    double d,
                    int n) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> x_rand(-d, d);
    std::uniform_real_distribution<double> y_rand(-2, 2);
//    std::uniform_real_distribution<double> d_rand(-d/10., d/10.);
    std::uniform_real_distribution<double> d_rand(0., 0.);

    n = std::max(n, 20);
    int id = 1;
    // side 1
    for (size_t i = 0; i < n; i++) {
        Eigen::Vector4d pt_tmp(x_rand(generator), x_rand(generator), 5 + d_rand(generator), id++);
        points.push_back(pt_tmp);
    }

    find3dPlanes_pts(points, plane_associate, plane_params);

    // save points
    save_points("all_points.txt", points);

}

void CreatePoints(Points &points) {
    std::ifstream f;
    f.open("house_model/house.txt");

    while (!f.eof()) {
        std::string s;
        std::getline(f, s);
        if (!s.empty()) {
            std::stringstream ss;
            ss << s;
            double x, y, z;
            ss >> x;
            ss >> y;
            ss >> z;
            Eigen::Vector4d pt0(x, y, z, 1);
            ss >> x;
            ss >> y;
            ss >> z;
            Eigen::Vector4d pt1(x, y, z, 1);

            pt0 = pt0 / 2.0;
            pt0(3) = 1.;
            pt1 = pt1 / 2.0;
            pt1(3) = 1.;

            bool isHistoryPoint = false;
            for (int i = 0; i < points.size(); ++i) {
                Eigen::Vector4d pt = points[i];
                if ((pt - pt0).norm() < 1e-3) {
                    isHistoryPoint = true;
                }
            }
            if (!isHistoryPoint)
                points.push_back(pt0);

            isHistoryPoint = false;
            for (int i = 0; i < points.size(); ++i) {
                Eigen::Vector4d pt = points[i];
                if ((pt - pt1).norm() < 1e-3) {
                    isHistoryPoint = true;
                }
            }
            if (!isHistoryPoint)
                points.push_back(pt1);
        }
    }

    // create more 3d points, you can comment this code
    int n = points.size();
    for (int j = 0; j < n; ++j) {
        Eigen::Vector4d p = points[j] + Eigen::Vector4d(0.5, 0.5, -0.5, 0);
        points.push_back(p);
    }

    // save points
    save_points("all_points.txt", points);
}

int main() {

    // 建立keyframe文件夹
    mkdir("keyframe", 0777);

    // default parameters
    Param params("../config.yaml");

    // 生成3d points
    Points points;
    std::map<int, int> plane_associate;
    std::vector<Eigen::Vector4d> plane_params;

    if (params.constraint_type == 0) {
        CreatePoints_wall(points, plane_associate, plane_params, 7, params.mpt_n);
    } else {
        CreatePoints_ground(points, plane_associate, plane_params, 7, params.mpt_n);
    }

    for (int i = 0; i < plane_params.size(); i++) {
        std::cout << "plane id: " << i << " parameter: " << plane_params[i].transpose() << std::endl;
    }

    for (const auto &id: plane_associate) {
        // std::cout << "pt id: " << id.first << " " << " plane id: " << id.second << std::endl;
    }


    /// generate cam pose
    std::vector<MotionData> camdata;
    std::cout << " start create cam data" << std::endl;
    {
        // step1. The first n camera poses remain still
        double static_time = 0.;
        int n = 10;
        for (int i = 0; i < n; i++) {
            MotionData cam;
            cam.timestamp = static_time;
            camdata.push_back(cam);
            static_time += 1.0 / params.cam_frequency;
        }

        // step2. camera moving                
        double moving_time = params.t_end;
        for (double t = 0; t < moving_time;) {

            // param
            double ellipse_x = 3;
            double ellipse_y = 4;
            double z = 0.5;                        // z-axis sin trajectory
            double K1 = 10;                      // z轴的正弦频率是 x，y 的 K1 倍
            double K = 2 * M_PI / moving_time;    // moving_time * K = 2pi 　　由于我们采取的是时间是20s, 系数K控制yaw正好旋转一圈，运动一周

            // translation
            // twb:  body frame in world frame
            Eigen::Vector3d tc0c(ellipse_x * cos(K * t), z * sin(K1 * K * t), ellipse_y * sin(K * t));
            // Rotation
            Eigen::Vector3d eulerAngles;
            if (params.constraint_type == 0) {
                eulerAngles = Eigen::Vector3d(0.1 * cos(t), -K * t + M_PI_2,
                                              0.2 * sin(t));   // roll ~ [-0.1, 0.1], pitch ~ [-0.2, 0.2], yaw ~ [0,2pi]
            } else {
                eulerAngles = Eigen::Vector3d(-1, -K * t + M_PI_2,
                                              0.2 * sin(t));   // roll ~ [-0.1, 0.1], pitch ~ [-0.2, 0.2], yaw ~ [0,2pi]
            }
            Eigen::Matrix3d Rc0c = euler2Rotation(eulerAngles);         // body frame to world frame



            Eigen::AngleAxisd rotYZ(90.0 / 180.0 * M_PI, Eigen::Vector3d::UnitX());

            MotionData cam;
            cam.timestamp = t + static_time;
            cam.Rwb = rotYZ.toRotationMatrix() * Rc0c;
            cam.twb = rotYZ.toRotationMatrix() * tc0c;

            camdata.push_back(cam);

            t += 1.0 / params.cam_frequency;
        }

        // step3. copy pose to the first n poses
        for (int i = 0; i < n; i++) {
            camdata[i].Rwb = camdata[n].Rwb;
            camdata[i].twb = camdata[n].twb;
        }

        save_Pose("cam_pose.txt", camdata);
        save_Pose_asTUM("cam_pose_tum.txt", camdata);
    }

    /// generate imu data with B-spline
    std::vector<MotionData> splinedata;
    std::cout << " start create spline data" << std::endl;
    {
        // IMU model
        IMU imuGen(params);

        imuGen.SplineMotionModel(camdata, splinedata, params.imu_frequency);

        if (params.add_imu_noise) {
            for (size_t i = 0; i < splinedata.size(); i++) {
                imuGen.addIMUnoise(splinedata[i]);
            }
        }

        save_Pose("spline_pose.txt", splinedata);
    }

    // points obs in image
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> noise(0.0, params.pixel_noise);

    // 按照SFM的数据结构来组织特征点
    std::unordered_map<FeatureID, SFMFeature> SFMConstruct;
    std::vector<Eigen::Matrix4d> Twcs;

    for (int n = 0; n < camdata.size(); ++n) {
        MotionData data = camdata[n];
        Eigen::Matrix4d Twc = Eigen::Matrix4d::Identity();
        Twc.block(0, 0, 3, 3) = data.Rwb;
        Twc.block(0, 3, 3, 1) = data.twb;
        Twcs.push_back(Twc);

        // 遍历所有的特征点，看哪些特征点在视野里
        std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d>> points_cam;    // ３维点在当前cam视野里
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> features_cam;  // 对应的２维图像坐标
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> features_normalize;      // 对应的２维图像坐标，用于 vins-mono
        for (int i = 0; i < points.size(); ++i) {
            Eigen::Vector4d pw = points[i];          // 最后一位存着feature id
            int id = int(pw[3]);
            pw[3] = 1;                               //改成齐次坐标最后一位
            Eigen::Vector4d pc1 = Twc.inverse() * pw; // T_wc.inverse() * Pw  -- > point in cam frame

            if (pc1(2) < 0) continue; // z必须大于０,在摄像机坐标系前方

            double u = pc1(0) / pc1(2) * params.fx + params.cx;
            double v = pc1(1) / pc1(2) * params.fy + params.cy;

            if (params.add_cam_noise) {
                u += noise(generator);
                v += noise(generator);
            }
            Eigen::Vector2d obs(u, v);
            if (u < params.image_w && u > 0 &&
                v > 0 && v < params.image_h) {
                points_cam.push_back(points[i]);
                features_cam.push_back(obs);

                Eigen::Vector2d obs_tmp((u - params.cx) / params.fx, (v - params.cy) / params.fy);
                features_normalize.push_back(obs_tmp);
            }
        }

        TimeFrameId frame_id = n;

        for (int i = 0; i < points_cam.size(); i++) {
            Eigen::Vector4d pt_world = points_cam[i];
            Eigen::Vector3d pt_normal = Eigen::Vector3d(features_normalize[i].x(), features_normalize[i].y(),
                                                        1.); // 归一化平面点
            Eigen::Vector2d pt_img = features_cam[i];    // 图像平面点

            int feature_id = (int) pt_world(3);
            int plane_id = plane_associate.at(feature_id);

            Eigen::Matrix4d Twc = Twcs[frame_id];
            Eigen::Vector4d pt_camera = Twc.inverse() * Eigen::Vector4d(pt_world.x(), pt_world.y(), pt_world.z(), 1);
            double Idepth_gt = 1 / pt_camera(2);

            FeaturePerFrame kpt_obs(pt_normal, pt_img);

            if (SFMConstruct.find(feature_id) == SFMConstruct.end()) {
                SFMConstruct[feature_id] = SFMFeature(feature_id, plane_id, frame_id);
                SFMConstruct[feature_id].Idepth_gt = Idepth_gt;
                SFMConstruct[feature_id].obs[frame_id] = kpt_obs;
            } else {
                SFMConstruct[feature_id].obs[frame_id] = kpt_obs;
            }
        }

        // save points
        std::stringstream filename1;
        filename1 << "keyframe/all_points_" << n << ".txt";
        save_features(filename1.str(), points_cam, features_cam);
        filename1.str(std::string());
        filename1.clear();
        filename1 << "keyframe/all_points_normal_" << n << ".txt";
        save_features(filename1.str(), points_cam, features_normalize);
    }

    std::cout << "number of points: " << SFMConstruct.size() << std::endl;
    std::cout << "number of poses: " << Twcs.size() << std::endl;


    BundleAdjustmentBase bundleAdjustment(params, SFMConstruct, Twcs, plane_params);

    // 1. triangulate
    bundleAdjustment.triangulate();
    // 2. add pose noise
    bundleAdjustment.add_pose_noise();
    // 3. add plane noise
    bundleAdjustment.add_wall_plane_noise();

    std::cout << "before optimization error: " << std::endl;
    evaluateTrajectory(Twcs, bundleAdjustment.Twcs);
    evaluatePlane(plane_params, bundleAdjustment.planes);


    TicToc t_solve_loose;
    bundleAdjustment.loosely_coupled_method();
    std::cout << "-------loose time used: " << t_solve_loose.toc() << " ms" << std::endl;

    std::cout << "after optimization error: " << std::endl;
    evaluateTrajectory(Twcs, bundleAdjustment.Twcs_opt);
    evaluatePlane(plane_params, bundleAdjustment.planes_opt);

    std::cout << "--------------------------------" << std::endl;

    TicToc t_solve_tight;
    bundleAdjustment.tightly_coupled_method();
    std::cout << "-------tight time used: " << t_solve_tight.toc() << " ms" << std::endl;

    std::cout << "after optimization error: " << std::endl;
    evaluateTrajectory(Twcs, bundleAdjustment.Twcs_opt);
    evaluatePlane(plane_params, bundleAdjustment.planes_opt);

    std::cout << "--------------------------------" << std::endl;

    TicToc t_solve_point;
    bundleAdjustment.point_method();
    std::cout << "-------point time used: " << t_solve_point.toc() << " ms" << std::endl;

    std::cout << "after optimization error: " << std::endl;
    evaluateTrajectory(Twcs, bundleAdjustment.Twcs_opt);

    viewer view(params);
    view.guiRun("cam_pose.txt", "all_points.txt", "keyframe/all_points_");
    return 0;
}
