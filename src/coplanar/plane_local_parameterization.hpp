//
// Created by xubo on 21-8-29.
//

#ifndef BASALT_LOCAL_PARAMETERIZATION_PLANE_H
#define BASALT_LOCAL_PARAMETERIZATION_PLANE_H

#include <math.h>
#include <ceres/ceres.h>

namespace mesh {

    class WallLocalParameterization : public ceres::LocalParameterization {
    public:
        virtual ~WallLocalParameterization() {}

        virtual bool Plus(double const *T_raw, double const *delta_raw,
                          double *T_plus_delta_raw) const {

            Eigen::Map<Eigen::Matrix<double, 4, 1> const> const plane(T_raw);
            Eigen::Map<Eigen::Vector3d const> const delta(delta_raw);
            Eigen::Map<Eigen::Matrix<double, 4, 1> > plane_plus_delta(T_plus_delta_raw);

            double d = plane.w();

            double delta_theta = delta(0);
            double delta_phi = delta(1);


            double theta= std::atan2(plane.y(), plane.x());
            double phi = std::asin(plane.z());
            theta += delta_theta;
            phi += delta_phi;
            Eigen::Vector3d n(cos(phi) * cos(theta), cos(phi) * sin(theta), sin(phi));

            d += delta(2);
            plane_plus_delta.head<3>() =  n;
            plane_plus_delta(3) = d;

            double norm = plane_plus_delta.head<3>().norm();
            plane_plus_delta = plane_plus_delta * (1. / norm);

            if (plane_plus_delta(3) > 0) {
                plane_plus_delta = -plane_plus_delta;
            }


            return true;
        }


        virtual bool ComputeJacobian(double const *T_raw,
                                     double *jacobian_raw) const {


            Eigen::Map<Eigen::Matrix<double, 4, 1> const> plane_w(T_raw);

            Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>> jacobian(
                    jacobian_raw);


            double theta = std::atan2(plane_w.y(), plane_w.x());
            double phi = std::asin(plane_w.z());
            jacobian.setZero();
            jacobian(0, 0) = -std::cos(phi) * std::sin(theta);
            jacobian(0, 1) = -std::sin(phi) * std::cos(theta);
            jacobian(1, 0) = std::cos(phi) * std::cos(theta);
            jacobian(1, 1) = -std::sin(phi) * std::sin(theta);
            jacobian(2, 1) = std::cos(phi);
            jacobian(3, 2) = 1.;

            return true;
        }

        virtual int GlobalSize() const { return 4; }

        virtual int LocalSize() const { return 3; }
    };

}
#endif //BASALT_LOCAL_PARAMETERIZATION_PLANE_H
