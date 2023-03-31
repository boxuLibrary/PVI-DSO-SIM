//
// Created by hyj on 17-6-22.
//

#include <iostream>
#include "param.h"
#include<opencv2/core/core.hpp>

Param::Param(std::string strSettingPath)
{
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);

    std::cout <<"--------- Load param: "<< strSettingPath <<" --------- \n"<< std::endl;

    cv::FileNode node = fSettings["IMU.Frequency"];
    imu_frequency = node.operator int();
    node = fSettings["Camera.fps"];
    cam_frequency = node.operator int();
    imu_timestep = 1./imu_frequency;
    cam_timestep = 1./cam_frequency;

    node = fSettings["Camera.fx"];
    fx = node.real();
    node = fSettings["Camera.fy"];
    fy = node.real();
    node = fSettings["Camera.cx"];
    cx = node.real();
    node = fSettings["Camera.cy"];
    cy = node.real();
    node = fSettings["Camera.width"];
    image_w = node.operator int();
    node = fSettings["Camera.height"];
    image_h = node.operator int();

    cv::Mat Tbc;
    node = fSettings["Tbc"];
    Tbc = node.mat();

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            R_bc(i,j) = Tbc.at<float>(i,j);
        }

    }
    t_bc = Eigen::Vector3d(Tbc.at<float>(0,3),Tbc.at<float>(1,3),Tbc.at<float>(2,3));
    // std::cout << R_bc << "\n" << t_bc << std::endl;

    node = fSettings["Camera.NoiseFlag"];
    add_cam_noise = node.operator int();
    node = fSettings["IMU.NoiseFlag"];
    add_imu_noise = node.operator int();

    node = fSettings["IMU.NoiseGyro"];
    gyro_noise_sigma = node.real();
    node = fSettings["IMU.NoiseAcc"];
    acc_noise_sigma = node.real();
    node = fSettings["IMU.GyroWalk"];
    gyro_bias_sigma = node.real();
    node = fSettings["IMU.AccWalk"];
    acc_bias_sigma = node.real();

    node = fSettings["IMU.GyroBiasInit"];
    gyro_bias_init = node.real();
    node = fSettings["IMU.AccBiasInit"];
    acc_bias_init = node.real();
    node = fSettings["Camera.PixelNoise"];
    pixel_noise = node.real();

    node = fSettings["Constraint.type"];
    constraint_type = node.real();

    node = fSettings["optimization_num"];
    optimization_num = node.real();

    std::cout << "IMU gyro noise: " << gyro_noise_sigma << " rad/s/sqrt(Hz)" << std::endl;
    std::cout << "IMU gyro walk: " << gyro_bias_sigma << " rad/s^2/sqrt(Hz)" << std::endl;
    std::cout << "IMU accelerometer noise: " << acc_noise_sigma << " m/s^2/sqrt(Hz)" << std::endl;
    std::cout << "IMU accelerometer walk: " << acc_bias_sigma << " m/s^3/sqrt(Hz)" << std::endl;
    std::cout << "IMU gyro bias init: " << gyro_bias_init << " m/s" << std::endl;
    std::cout << "IMU accelerometer bias init: " << acc_bias_init << " m/s^2" << std::endl;

    node = fSettings["nFeatures"];
    mpt_n = node.operator int();
    node = fSettings["Times"];
    t_end = node.real();

    std::cout << "Map points number: " << mpt_n << std::endl;
    std::cout << "add cam noise flag: " << add_cam_noise << " imu noise: " << add_imu_noise << std::endl;

    std::cout << "add constraint type: " << constraint_type << std::endl;

    std::cout <<"\n--------- Load param END --------- "<< std::endl;
}
