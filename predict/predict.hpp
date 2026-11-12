#pragma once

#include <tuple>

#include <Eigen/Dense>

#include "ean_mean.hpp"
#include "ecsrt_v_ist.hpp"
#include "ist_v_ecrt.hpp"
#include "kep_dec.hpp"
#include "mean_motion.hpp"
#include "solve_kepler.hpp"

using State = std::tuple<Eigen::Vector3d, Eigen::Vector3d>;

State predict(const Eigen::Vector3d &r0, const Eigen::Vector3d &v0, double t0, double t, double mu);