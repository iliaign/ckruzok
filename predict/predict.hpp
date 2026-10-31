// predict.hpp - мастер-заголовок библиотеки predict_lib
#pragma once

#include <tuple>

#include <Eigen/Dense>

#include "kep_dec.hpp"
#include "ist_v_ecrt.hpp"
#include "ean_mean.hpp"
#include "mean_motion.hpp"
#include "solve_kepler.hpp"
#include "ecsrt_v_ist.hpp"

std::tuple<Eigen::Vector3d, Eigen::Vector3d> predict(
    const Eigen::Vector3d& r0,
    const Eigen::Vector3d& v0,
    double t0,
    double t,
    double mu
    );