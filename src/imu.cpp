//
// Created by hyj on 18-1-19.
//

#include <random>
#include "imu.h"
#include "utilities.h"
#include "se3.hpp"
#include "so3.hpp"

// euler2Rotation:   body frame to inertial frame
Eigen::Matrix3d euler2Rotation( Eigen::Vector3d  eulerAngles)
{
    double roll = eulerAngles(0);
    double pitch = eulerAngles(1);
    double yaw = eulerAngles(2);

    double cr = cos(roll); double sr = sin(roll);
    double cp = cos(pitch); double sp = sin(pitch);
    double cy = cos(yaw); double sy = sin(yaw);

    Eigen::Matrix3d RIb;
    RIb<< cy*cp ,   cy*sp*sr - sy*cr,   sy*sr + cy* cr*sp,
            sy*cp,    cy *cr + sy*sr*sp,  sp*sy*cr - cy*sr,
            -sp,         cp*sr,           cp*cr;
    return RIb;
}

Eigen::Matrix3d eulerRates2bodyRates(Eigen::Vector3d eulerAngles)
{
    double roll = eulerAngles(0);
    double pitch = eulerAngles(1);

    double cr = cos(roll); double sr = sin(roll);
    double cp = cos(pitch); double sp = sin(pitch);

    Eigen::Matrix3d R;
    R<<  1,   0,    -sp,
            0,   cr,   sr*cp,
            0,   -sr,  cr*cp;

    return R;
}


IMU::IMU(Param p): param_(p)
{
    gyro_bias_ = Eigen::Vector3d(1, 1, 1) * param_.gyro_bias_init;
    acc_bias_ = Eigen::Vector3d(1,1,1) * param_.acc_bias_init;
}

void IMU::addIMUnoise(MotionData& data)
{
//    std::random_device rd;
//    std::default_random_engine generator_(rd());
//    std::normal_distribution<double> noise(0.0, 1.0);
//
//    Eigen::Vector3d noise_gyro(noise(generator_),noise(generator_),noise(generator_));
//    Eigen::Matrix3d gyro_sqrt_cov = param_.gyro_noise_sigma * Eigen::Matrix3d::Identity();
//    data.imu_gyro = data.imu_gyro + gyro_sqrt_cov * noise_gyro / sqrt( param_.imu_timestep ) + gyro_bias_;
//
//    Eigen::Vector3d noise_acc(noise(generator_),noise(generator_),noise(generator_));
//    Eigen::Matrix3d acc_sqrt_cov = param_.acc_noise_sigma * Eigen::Matrix3d::Identity();
//    data.imu_acc = data.imu_acc + acc_sqrt_cov * noise_acc / sqrt( param_.imu_timestep ) + acc_bias_;
//
//    // gyro_bias update
//    Eigen::Vector3d noise_gyro_bias(noise(generator_),noise(generator_),noise(generator_));
//    gyro_bias_ += param_.gyro_bias_sigma * sqrt(param_.imu_timestep ) * noise_gyro_bias;
//    data.imu_gyro_bias = gyro_bias_;
//
//    // acc_bias update
//    Eigen::Vector3d noise_acc_bias(noise(generator_),noise(generator_),noise(generator_));
//    acc_bias_ += param_.acc_bias_sigma * sqrt(param_.imu_timestep ) * noise_acc_bias;
//    data.imu_acc_bias = acc_bias_;

    gyro_bias_ = Eigen::Vector3d(0.1, 0.1, 0.1);
    data.imu_gyro = data.imu_gyro + gyro_bias_;

}

MotionData IMU::MotionModel(double t)
{

    MotionData data;
    // param
    float ellipse_x = 7.5;
    float ellipse_y = 10;
    float z = 1;           // z轴做sin运动
    float K1 = 10;          // z轴的正弦频率是x，y的k1倍
    float K = M_PI/ 10;    // 20 * K = 2pi 　　由于我们采取的是时间是20s, 系数K控制yaw正好旋转一圈，运动一周

    // translation
    // twb:  body frame in world frame
    Eigen::Vector3d position( ellipse_x * cos( K * t) + 2.5, ellipse_y * sin( K * t) + 2.5,  z * sin( K1 * K * t ) + 2.5);
    Eigen::Vector3d dp(- K * ellipse_x * sin(K*t),  K * ellipse_y * cos(K*t), z*K1*K * cos(K1 * K * t));              // position导数　in world frame
    double K2 = K*K;
    Eigen::Vector3d ddp( -K2 * ellipse_x * cos(K*t),  -K2 * ellipse_y * sin(K*t), -z*K1*K1*K2 * sin(K1 * K * t));     // position二阶导数

    // Rotation
    double k_roll = 0.1;
    double k_pitch = 0.2;
    Eigen::Vector3d eulerAngles(k_roll * cos(t) , k_pitch * sin(t) , K*t );   // roll ~ [-0.2, 0.2], pitch ~ [-0.3, 0.3], yaw ~ [0,2pi]
    Eigen::Vector3d eulerAnglesRates(-k_roll * sin(t) , k_pitch * cos(t) , K);      // euler angles 的导数

//    Eigen::Vector3d eulerAngles(0.0,0.0, K*t );   // roll ~ 0, pitch ~ 0, yaw ~ [0,2pi]
//    Eigen::Vector3d eulerAnglesRates(0.,0. , K);      // euler angles 的导数

    Eigen::Matrix3d Rwb = euler2Rotation(eulerAngles);         // body frame to world frame
    Eigen::Vector3d imu_gyro = eulerRates2bodyRates(eulerAngles) * eulerAnglesRates;   //  euler rates trans to body gyro

    Eigen::Vector3d gn (0,0,-9.81);                                   //  gravity in navigation frame(ENU)   ENU (0,0,-9.81)  NED(0,0,9,81)
    Eigen::Vector3d imu_acc = Rwb.transpose() * ( ddp -  gn );  //  Rbw * Rwn * gn = gs

    data.imu_gyro = imu_gyro;
    data.imu_acc = imu_acc;
    data.Rwb = Rwb;
    data.twb = position;
    data.imu_velocity = dp;
    data.timestamp = t;
    return data;

}

/* Alonso Patron-Perez, 2015, A Spline-Based Trajectory Representation for Sensor Fusion and Rolling Shutter Cameras.
 * Patrick Geneva, 2018, https://udel.edu/~pgeneva/downloads/notes/2018_notes_mueffler2017arxiv.pdf
*/
void IMU::SplineMotionModel(std::vector< MotionData > cam_traj, std::vector< MotionData >& spline_imus, int frequency)
{
    if(cam_traj.size() < 4) return;

    Eigen::Vector3d gn (0,0,-9.81);            //  gravity in navigation frame(ENU)   ENU (0,0,-9.81)  NED(0,0,9,81)
    double start_time = cam_traj[1].timestamp;
    int cam_n = cam_traj.size() - 2;
    double end_time = cam_traj[cam_n].timestamp;
    double dt = 1.0/ frequency;
    // create imu data with timestamp
    for (double t = start_time; t < end_time + 1e-5; t += dt)
    {
        MotionData imu_tmp;
        imu_tmp.timestamp = t;
        spline_imus.push_back(imu_tmp);
    }
    // std::cout << cam_n << std::endl;

    int imu_index = 0;
    for (size_t i = 1; i < cam_traj.size()-2; i++)
    {
        double ti = cam_traj[i].timestamp;
        double tj = cam_traj[i+1].timestamp;

        // 4 control pose: i-1, i, i+1, i+2
        std::vector< Sophus::SE3d > Twbs;
        for (int j = -1; j < 3; j++)
        {
            Eigen::Matrix3d Rwb = cam_traj[i+j].Rwb * param_.R_bc.transpose();
            Eigen::Vector3d twb = -Rwb * param_.t_bc + cam_traj[i+j].twb;
            Sophus::SE3d Twb(Rwb, twb);
            Twbs.push_back(Twb);
        }

        Sophus::Vector6d Omege_im1_i   = (Twbs[0].inverse() * Twbs[1]).log();
        Sophus::Vector6d Omege_i_ip1   = (Twbs[1].inverse() * Twbs[2]).log();
        Sophus::Vector6d Omege_ip1_ip2 = (Twbs[2].inverse() * Twbs[3]).log();

        double delta_t = tj - ti;

        // sample imu data from ti to tj
        while( imu_index < spline_imus.size() && spline_imus.at(imu_index).timestamp < tj)
        {
            double timu = spline_imus.at(imu_index).timestamp;
            if( timu >= ti && timu < tj)
            {
                // sample imu data from spline
                double u = (timu - ti)/delta_t;

                double B_0 = 1.0 / 6.0 * (5 + 3*u - 3*u*u + 1*u*u*u);
                double B_1 = 1.0 / 6.0 * (1 + 3*u + 3*u*u - 2*u*u*u);
                double B_2 = 1.0 / 6.0 * (u*u*u);

                double dB_0 = 1.0 / (6 * delta_t) * (3 - 6*u + 3*u*u);
                double dB_1 = 1.0 / (6 * delta_t) * (3 + 6*u - 6*u*u);
                double dB_2 = 1.0 / (6 * delta_t) * (3*u*u);

                double ddB_0 = 1.0 / (6 * delta_t * delta_t) * (-6 + 6*u);
                double ddB_1 = 1.0 / (6 * delta_t * delta_t) * (6 - 12*u);
                double ddB_2 = 1.0 / (6 * delta_t * delta_t) * (6*u);

                Eigen::Matrix<double, 4, 4> A_0 = Sophus::SE3d::exp(B_0 * Omege_im1_i).matrix();
                Eigen::Matrix<double, 4, 4> A_1 = Sophus::SE3d::exp(B_1 * Omege_i_ip1).matrix();
                Eigen::Matrix<double, 4, 4> A_2 = Sophus::SE3d::exp(B_2 * Omege_ip1_ip2).matrix();

                Eigen::Matrix<double, 4, 4> dA_0 = dB_0 * Sophus::SE3d::hat(Omege_im1_i).matrix() * A_0;
                Eigen::Matrix<double, 4, 4> dA_1 = dB_1 * Sophus::SE3d::hat(Omege_i_ip1).matrix() * A_1;
                Eigen::Matrix<double, 4, 4> dA_2 = dB_2 * Sophus::SE3d::hat(Omege_ip1_ip2).matrix() * A_2;
                // Eigen::Matrix<double, 4, 4> dA_1 = dB_1 * A_1 * Sophus::SE3::hat(Omege_i_ip1).matrix();
                // Eigen::Matrix<double, 4, 4> dA_2 = dB_2 * A_2 * Sophus::SE3::hat(Omege_ip1_ip2).matrix();

                Eigen::Matrix<double, 4, 4> ddA_0 = dB_0 * Sophus::SE3d::hat(Omege_im1_i).matrix() * dA_0 + ddB_0 * Sophus::SE3d::hat(Omege_im1_i).matrix() * A_0;
                Eigen::Matrix<double, 4, 4> ddA_1 = dB_1 * Sophus::SE3d::hat(Omege_i_ip1).matrix() * dA_1 + ddB_1 * Sophus::SE3d::hat(Omege_i_ip1).matrix() * A_1;
                Eigen::Matrix<double, 4, 4> ddA_2 = dB_2 * Sophus::SE3d::hat(Omege_ip1_ip2).matrix() * dA_2 + ddB_2 * Sophus::SE3d::hat(Omege_ip1_ip2).matrix() * A_2;
                // Eigen::Matrix<double, 4, 4> ddA_1 = dB_1 * dA_1 * Sophus::SE3::hat(Omege_i_ip1).matrix() + ddB_1 * A_1 * Sophus::SE3::hat(Omege_i_ip1).matrix();
                // Eigen::Matrix<double, 4, 4> ddA_2 = dB_2 * dA_2 * Sophus::SE3::hat(Omege_ip1_ip2).matrix() + ddB_2 * A_2 * Sophus::SE3::hat(Omege_ip1_ip2).matrix();

                Eigen::Matrix<double, 4, 4> T_w_s = Twbs[0].matrix() * A_0 * A_1 * A_2;
                Eigen::Matrix<double, 4, 4> dT_w_s = Twbs[0].matrix() * (dA_0*A_1*A_2 + A_0*dA_1*A_2 + A_0*A_1*dA_2);
                Eigen::Matrix<double, 4, 4> ddT_w_s = Twbs[0].matrix() * (ddA_0*A_1*A_2 + A_0*ddA_1*A_2 + A_0*A_1*ddA_2 + 2*(dA_0*dA_1*A_2 + A_0*dA_1*dA_2 + dA_0*A_1*dA_2));

                Eigen::Matrix<double, 3, 3> R_w_s = T_w_s.block(0,0,3,3);
                Eigen::Matrix<double, 3, 3> dR_w_s = dT_w_s.block(0,0,3,3);
                Eigen::Vector3d dp_w_s = dT_w_s.block<3,1>(0,3);
                Eigen::Vector3d ddp_w_s = ddT_w_s.block<3,1>(0,3);

                spline_imus.at(imu_index).Rwb = R_w_s;
                spline_imus.at(imu_index).twb = T_w_s.block<3,1>(0,3);
                spline_imus.at(imu_index).imu_gyro = Sophus::SO3d::vee(R_w_s.transpose() * dR_w_s);  // dR = R * [omega_b]_{\times} --> [omega_b]_{\times} = R^{\top} * dR
                spline_imus.at(imu_index).imu_acc = R_w_s.transpose() * (ddp_w_s - gn);
                spline_imus.at(imu_index).imu_velocity = dp_w_s;
            }

            imu_index++;
        }

    }

}

//读取生成的imu数据并用imu动力学模型对数据进行计算，最后保存imu积分以后的轨迹，
//用来验证数据以及模型的有效性。
void IMU::testImu(std::string src, std::string dist)
{
    std::vector<MotionData>imudata;
    LoadPose(src,imudata);

    std::ofstream save_points;
    save_points.open(dist);

    double dt = param_.imu_timestep;
    Eigen::Vector3d Pwb = init_twb_;              // position :    from  imu measurements
    Eigen::Quaterniond Qwb(init_Rwb_);            // quaterniond:  from imu measurements
    Eigen::Vector3d Vw = init_velocity_;          // velocity  :   from imu measurements
    Eigen::Vector3d gw(0,0,-9.81);    // ENU frame
    Eigen::Vector3d temp_a;
    Eigen::Vector3d theta;
    for (int i = 1; i < imudata.size(); ++i) {

        MotionData imupose = imudata[i];

        //delta_q = [1 , 1/2 * thetax , 1/2 * theta_y, 1/2 * theta_z]
        Eigen::Quaterniond dq;
        Eigen::Vector3d dtheta_half =  imupose.imu_gyro * dt /2.0;
        dq.w() = 1;
        dq.x() = dtheta_half.x();
        dq.y() = dtheta_half.y();
        dq.z() = dtheta_half.z();
        dq.normalize();

        /// imu 动力学模型 欧拉积分
        Eigen::Vector3d acc_w = Qwb * (imupose.imu_acc) + gw;  // aw = Rwb * ( acc_body - acc_bias ) + gw
        Qwb = Qwb * dq;
        Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
        Vw = Vw + acc_w * dt;

        /// 中值积分

        //　按着imu postion, imu quaternion , cam postion, cam quaternion 的格式存储，由于没有cam，所以imu存了两次
        save_points<<imupose.timestamp<<" "
                   <<Qwb.w()<<" "
                   <<Qwb.x()<<" "
                   <<Qwb.y()<<" "
                   <<Qwb.z()<<" "
                   <<Pwb(0)<<" "
                   <<Pwb(1)<<" "
                   <<Pwb(2)<<" "
                   <<Qwb.w()<<" "
                   <<Qwb.x()<<" "
                   <<Qwb.y()<<" "
                   <<Qwb.z()<<" "
                   <<Pwb(0)<<" "
                   <<Pwb(1)<<" "
                   <<Pwb(2)<<" "
                   <<std::endl;

    }

    std::cout<<"test　end"<<std::endl;

}
