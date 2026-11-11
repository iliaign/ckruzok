// kep_to_dec.hpp
#pragma once

#include <tuple>

#include "dec_to_kep.hpp"
#include <Eigen/Dense>

std::tuple<Eigen::Vector3d, Eigen::Vector3d> convert(const Kepler &kep, double mu);
