// dec_to_kep.cpp
#include "dec_to_kep.hpp"
#include <cmath>

static double wrap_angle(const double x_in) {
    double x = std::fmod(x_in, 2.0 * M_PI);
    if (x < 0.0) x += 2.0 * M_PI;
    return x;
}

static double signed_angle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& n) {
    return std::atan2(a.cross(b).dot(n), a.dot(b));
}

Kepler convert(const Eigen::Vector3d& r,const Eigen::Vector3d& v,  const double mu) {
    const double rn = r.norm();
    const double vn2 = v.squaredNorm();
    const Eigen::Vector3d l = r.cross(v);
    const Eigen::Vector3d l_hat = l.normalized();
    const Eigen::Vector3d n_raw = Eigen::Vector3d::UnitZ().cross(l);
    const bool equatorial = (n_raw.squaredNorm() == 0.0);
    const Eigen::Vector3d nvec = equatorial ? Eigen::Vector3d(Eigen::Vector3d::UnitX()) : n_raw;
    const Eigen::Vector3d evec = ((vn2 - mu / rn) * r - r.dot(v) * v) / mu;
    //вторая часть скобки, если ортиа круговая то вектор
    //эксентриситета совпадает с напрвелнеим на 
    //восходящий узел
    const bool circular = (evec.squaredNorm() == 0.0);
    const Eigen::Vector3d e1 = circular ? nvec : evec;

    const double inv_a = 2.0 / rn - vn2 / mu;
    const double a =  1.0 / inv_a;

    return Kepler{
        a,                                                    // большая полуось
        evec.norm(),                                          // эксцентриситет
        std::atan2(l.head<2>().norm(), l.z()),                // наклонение орбиты
        wrap_angle(std::atan2(nvec.y(), nvec.x())),           // долгота восходящего узла
        wrap_angle(signed_angle(nvec, e1, l_hat)),            // аргумент перицентра
        wrap_angle(signed_angle(e1, r, l_hat))                // истинная аномалия
    };
}

