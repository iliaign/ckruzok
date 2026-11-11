// dec_to_kep.hpp
#pragma once

#include "kepler.hpp"
#include <Eigen/Dense>

Kepler convert(const Eigen::Vector3d &r, const Eigen::Vector3d &v, double mu);