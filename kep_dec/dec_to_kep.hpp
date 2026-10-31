// dec_to_kep.hpp
#pragma once

#include <Eigen/Dense>
#include "kepler.hpp"

Kepler convert(const Eigen::Vector3d& r, const Eigen::Vector3d& v, double mu);