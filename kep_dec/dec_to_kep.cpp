// dec_to_kep.cpp
#include "dec_to_kep.hpp"

#include <algorithm>
#include <cmath>

// приведение угла к [0, 2pi)
static double wrap_angle(double x) {
    x = std::fmod(x, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x;
}

// знаковый угол от a к b вокруг единичной оси n (через atan2 - без NaN)
static double signed_angle(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                           const Eigen::Vector3d& n) {
    return std::atan2(a.cross(b).dot(n), a.dot(b));
}

Eigen::Matrix<double, 6, 1> dec_to_kep(const Eigen::Vector3d& r,
                                       const Eigen::Vector3d& v,
                                       double mu) {
    const double rn = r.norm();
    const double vn = v.norm();

    const Eigen::Vector3d h = r.cross(v);                              // момент импульса
    const Eigen::Vector3d h_hat = h.normalized();
    const Eigen::Vector3d nvec = Eigen::Vector3d::UnitZ().cross(h);    // линия узлов
    const Eigen::Vector3d evec = ((vn * vn - mu / rn) * r - r.dot(v) * v) / mu;

    Eigen::Matrix<double, 6, 1> kep;
    kep(0) = 1.0 / (2.0 / rn - vn * vn / mu);                          // a
    kep(1) = evec.norm();                                              // e
    kep(2) = std::acos(std::clamp(h_hat.z(), -1.0, 1.0));              // i
    kep(3) = wrap_angle(std::atan2(nvec.y(), nvec.x()));               // raan
    kep(4) = wrap_angle(signed_angle(nvec, evec, h_hat));              // argp
    kep(5) = wrap_angle(signed_angle(evec, r, h_hat));                 // nu
    return kep;
}
