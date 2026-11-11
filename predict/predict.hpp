// predict.hpp
#pragma once

#include <Eigen/Dense>
#include <tuple>


using State = std::tuple<Eigen::Vector3d, Eigen::Vector3d>;

State predict(const Eigen::Vector3d &r0, const Eigen::Vector3d &v0, double t0, double t, double mu);