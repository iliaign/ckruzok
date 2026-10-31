// dec_to_kep.hpp
#pragma once

#include <Eigen/Dense>

namespace kep_dec {


struct KeplerOut {
    double a{};
    double e{};
    double inc{};
    double raan{};
    double argp{};
    double nu{};
};

// Декартовы r, v -> элементы Кеплера.
// Если mu == 0 или |r| == 0, возвращается KeplerOut{} (все нули).
KeplerOut convert(const Eigen::Vector3d& r,
                  const Eigen::Vector3d& v,
                  double mu);

}  // namespace kep_dec
