// dec_to_kep.cpp
#include "dec_to_kep.hpp"
#include <cmath>

// приведение угла к [0, 2pi)
double DecToKep::wrap_angle(double x) {
    x = std::fmod(x, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x;
}

double DecToKep::signed_angle(const Eigen::Vector3d& a,
                              const Eigen::Vector3d& b,
                              const Eigen::Vector3d& n) {
    return std::atan2(a.cross(b).dot(n), a.dot(b));
}

Eigen::Matrix<double, 6, 1> DecToKep::convert(const Eigen::Vector3d& r,
                                              const Eigen::Vector3d& v) const {

    const double rn = r.norm();
    const double vn = v.norm();

    //вектор орбиального момента 
    const Eigen::Vector3d l = r.cross(v);                              // орбитальный момент
    const Eigen::Vector3d l_hat = l.normalized();

    // экваториальная орбита линия узлов не определен
    // направление на восходящий узел совпадает с осью x
    Eigen::Vector3d nvec = Eigen::Vector3d::UnitZ().cross(l);          // линия узлов
    const bool equatorial = nvec.norm() < 1e-10;
    nvec = equatorial ? Eigen::Vector3d(Eigen::Vector3d::UnitX()) : nvec;

    const Eigen::Vector3d evec = ((vn * vn - mu_ / rn) * r - r.dot(v) * v) / mu_; //векторный эксентриситет
    // круговая орбита e1 = e/|e| при e != 0, иначе e1 = N 
    const Eigen::Vector3d e1 = (evec.norm() > 1e-10) ? evec.normalized() : nvec.normalized();

    Eigen::Matrix<double, 6, 1> kep;
    kep(0) = 1.0 / (2.0 / rn - vn * vn / mu_);                         // a

    kep(1) = evec.norm();                                              // e

    kep(2) = std::atan2( l.head<2>().norm(), l.z() );                  // L

    kep(3) = wrap_angle(std::atan2(nvec.y(), nvec.x()));               // долгтоа восходящего узла

    kep(4) = wrap_angle(signed_angle(nvec, e1, l_hat));                // аргумент перицентра 

    kep(5) = wrap_angle(signed_angle(e1, r, l_hat));                   // истинная анамалия

    return kep;
}