// dec_to_kep.hpp
#pragma once

#include <Eigen/Dense>

class DecToKep {
public:
    explicit DecToKep(double mu) : mu_(mu) {}

    Eigen::Matrix<double, 6, 1> convert(const Eigen::Vector3d& r,
                                        const Eigen::Vector3d& v) const;

private:
    static double wrap_angle(double x);
    static double signed_angle(const Eigen::Vector3d& a,
                               const Eigen::Vector3d& b,
                               const Eigen::Vector3d& n);

    double mu_;
};