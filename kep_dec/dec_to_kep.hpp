// dec_to_kep.hpp
// Декартово состояние (r, v) -> кеплеровские элементы.
// Элементы - вектор Eigen из 6 компонент (углы в радианах):
//   kep(0) = a     большая полуось
//   kep(1) = e     эксцентриситет
//   kep(2) = i     наклонение
//   kep(3) = raan  долгота восходящего узла
//   kep(4) = argp  аргумент перицентра
//   kep(5) = nu    истинная аномалия
#pragma once

#include <Eigen/Dense>

Eigen::Matrix<double, 6, 1> dec_to_kep(const Eigen::Vector3d& r,
                                       const Eigen::Vector3d& v,
                                       double mu);
