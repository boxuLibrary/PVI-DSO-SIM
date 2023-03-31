//
// Created by hyj on 17-6-22.
//

#ifndef IMUSIM_PARAM_H
#define IMUSIM_PARAM_H

#include <Eigen/Core>

class Param{

public:

    Param(std::string strSettingPath);

    // time
    int imu_frequency = 200;
    int cam_frequency = 10;
    double imu_timestep = 1./imu_frequency;
    double cam_timestep = 1./cam_frequency;
    double t_start = 0.;
    double t_end = 40;      // senconds

    // map point number
    int mpt_n = 100;

    // noise
    double gyro_bias_sigma = 1.0e-5;     // ad/s^2/sqrt(Hz)
    double acc_bias_sigma = 0.0001;      // m/s^3/sqrt(Hz)

    double gyro_noise_sigma = 0.015;     // rad/s * 1/sqrt(hz)
    double acc_noise_sigma = 0.019;      //　m/(s^2) * 1/sqrt(hz)

    double gyro_bias_init = 0.001;
    double acc_bias_init = 0.01;

    double pixel_noise = 0.1;            // 1 pixel noise

    int add_cam_noise = 0;
    int add_imu_noise = 0;

    // cam f
    double fx = 460;
    double fy = 460;
    double cx = 320;
    double cy = 240;
    int image_w = 640;
    int image_h = 480;

    int constraint_type = 0;

    int optimization_num = 10;


    // 外参数
    Eigen::Matrix3d R_bc;   // cam to body
    Eigen::Vector3d t_bc;     // cam to body

};


#endif //IMUSIM_PARAM_H