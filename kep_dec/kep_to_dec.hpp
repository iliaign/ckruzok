// kep_to_dec.hpp
#pragma once

#include <Eigen/Dense>

class KepToDec {
public:
    explicit KepToDec(double mu) : mu_(mu) {}

    void convert(const Eigen::Matrix<double, 6, 1>& kep,
                 Eigen::Vector3d& r, Eigen::Vector3d& v) const;

private:
    double mu_;
};