// kep_to_dec.hpp
#pragma once

#include <tuple>

#include <Eigen/Dense>

namespace kep_dec {

// Вход kep_to_dec: элементы Кеплера (углы в радианах)
struct KeplerIn {
    double a{};
    double e{};
    double inc{};
    double raan{};
    double argp{};
    double nu{};
};

// Элементы Кеплера -> (r, v).
// Если p == 0 или 1 + e*cos(nu) == 0, возвращаются нулевые векторы.
std::tuple<Eigen::Vector3d, Eigen::Vector3d> convert(const KeplerIn& kep,
                                                     double mu);

}  // namespace kep_dec
