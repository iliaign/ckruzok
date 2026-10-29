#include "kep_to_dec.hpp"

#include <cmath>

void KepToDec::convert(const Eigen::Matrix<double, 6, 1>& kep,
                       Eigen::Vector3d& r, Eigen::Vector3d& v) const {
    const double a    = kep(0);
    const double e    = kep(1);
    const double inc  = kep(2);
    const double raan = kep(3);
    const double argp = kep(4);
    const double nu   = kep(5);

    const double p   = a * (1.0 - e * e);       
    const double cnu = std::cos(nu);
    const double snu = std::sin(nu);
    const double rn  = p / (1.0 + e * cnu);
    const double sq  = std::sqrt(mu_ / p);

    const Eigen::Vector3d r_pqw(rn * cnu, rn * snu, 0.0);
    const Eigen::Vector3d v_pqw(-sq * snu, sq * (e + cnu), 0.0);

    const Eigen::Matrix3d Q =
        (Eigen::AngleAxisd(raan, Eigen::Vector3d::UnitZ()) *
         Eigen::AngleAxisd(inc,  Eigen::Vector3d::UnitX()) *
         Eigen::AngleAxisd(argp, Eigen::Vector3d::UnitZ())).toRotationMatrix();

    r = Q * r_pqw;
    v = Q * v_pqw;
}