// dec_to_kep.cpp
#include "dec_to_kep.hpp"

#include <cmath>
#include <limits>

namespace kep_dec {

namespace {

// приведение угла к [0, 2pi)
double wrap_angle(const double x_in) {
    double x = std::fmod(x_in, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x;
}

// n должен быть единичным (или нулевым); длины a и b на atan2 не влияют
double signed_angle(const Eigen::Vector3d& a,
                    const Eigen::Vector3d& b,
                    const Eigen::Vector3d& n) {
    return std::atan2(a.cross(b).dot(n), a.dot(b));
}

// нормировка без деления на ноль: нулевой вектор остаётся нулевым
Eigen::Vector3d unit_or_zero(const Eigen::Vector3d& x) {
    const double n = x.norm();
    return (n == 0.0) ? Eigen::Vector3d(Eigen::Vector3d::Zero())
                      : Eigen::Vector3d(x / n);
}

}  // namespace

KeplerOut convert(const Eigen::Vector3d& r,
               const Eigen::Vector3d& v,
               const double mu) {
    const double rn = r.norm();

    // делим на mu и на rn
    if (mu == 0.0 || rn == 0.0) {
        return KeplerOut{};
    }

    const double vn2 = v.squaredNorm();

    // вектор орбитального момента (нулевой для прямолинейной орбиты)
    const Eigen::Vector3d l = r.cross(v);
    const Eigen::Vector3d l_hat = unit_or_zero(l);

    // линия узлов; нулевая -> экваториальная орбита, берём ось X
    const Eigen::Vector3d n_raw = Eigen::Vector3d::UnitZ().cross(l);
    const bool equatorial = (n_raw.squaredNorm() == 0.0);
    const Eigen::Vector3d nvec =
        equatorial ? Eigen::Vector3d(Eigen::Vector3d::UnitX()) : n_raw;

    // вектор эксцентриситета (деление на mu != 0)
    const Eigen::Vector3d evec = ((vn2 - mu / rn) * r - r.dot(v) * v) / mu;

    // нулевой -> круговая орбита, направление отсчёта = узел
    const bool circular = (evec.squaredNorm() == 0.0);
    const Eigen::Vector3d e1 = circular ? nvec : evec;

    // 1/a == 0 -> параболическая орбита, a = inf
    const double inv_a = 2.0 / rn - vn2 / mu;
    const double a = (inv_a == 0.0) ? std::numeric_limits<double>::infinity()
                                    : 1.0 / inv_a;

    return KeplerOut{
        a,
        evec.norm(),                                        // e
        std::atan2(l.head<2>().norm(), l.z()),              // inc
        wrap_angle(std::atan2(nvec.y(), nvec.x())),         // raan
        wrap_angle(signed_angle(nvec, e1, l_hat)),          // argp
        wrap_angle(signed_angle(e1, r, l_hat))              // nu
    };
}

}  // namespace kep_dec
