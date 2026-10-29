// kep_to_dec.hpp
// Кеплеровские элементы -> декартово состояние (r, v).
// Порядок элементов: [a, e, i, raan, argp, nu], углы в радианах.
#pragma once

#include <Eigen/Dense>

void kep_to_dec(const Eigen::Matrix<double, 6, 1>& kep, double mu,
                Eigen::Vector3d& r, Eigen::Vector3d& v);
