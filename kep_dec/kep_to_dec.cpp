// kep_to_dec.cpp
#include "kep_to_dec.hpp"

#include <cmath>

namespace kep_dec {

std::tuple<Eigen::Vector3d, Eigen::Vector3d> convert(const KeplerIn& kep,
                                                     const double mu) {
    const double p     = kep.a * (1.0 - kep.e * kep.e);
    const double cnu   = std::cos(kep.nu);
    const double snu   = std::sin(kep.nu);
    const double denom = 1.0 + kep.e * cnu;

    // делим на denom (p / denom) и на p (mu / p)
    if (p == 0.0 || denom == 0.0) {
        return std::tuple<Eigen::Vector3d, Eigen::Vector3d>(
            Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    }

    const double rn = p / denom;
    const double sq = std::sqrt(mu / p);

    const Eigen::Vector3d r_pqw(rn * cnu, rn * snu, 0.0);
    const Eigen::Vector3d v_pqw(-sq * snu, sq * (kep.e + cnu), 0.0);

    const Eigen::Matrix3d Q =
        (Eigen::AngleAxisd(kep.raan, Eigen::Vector3d::UnitZ()) *
         Eigen::AngleAxisd(kep.inc,  Eigen::Vector3d::UnitX()) *
         Eigen::AngleAxisd(kep.argp, Eigen::Vector3d::UnitZ())).toRotationMatrix();

    return std::tuple<Eigen::Vector3d, Eigen::Vector3d>(Q * r_pqw, Q * v_pqw);
}

}  // namespace kep_dec
