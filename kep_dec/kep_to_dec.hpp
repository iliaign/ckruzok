// kep_to_dec.hpp
#pragma once

#include <tuple>

#include <Eigen/Dense>
#include "dec_to_kep.hpp"

std::tuple<Eigen::Vector3d, Eigen::Vector3d> convert(const KeplerIn& kep, double mu);

