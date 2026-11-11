// predict.cpp
#include "predict.hpp"

State predict(const Eigen::Vector3d &r0, const Eigen::Vector3d &v0, const double t0, const double t, const double mu)
{
    const Kepler k0 = convert(r0, v0, mu);
    const double E0 = ist_v_ecrt(k0.nu, k0.e);
    const double M0 = ean_mean(E0, k0.e);
    const double n = mean_motion(k0.a, mu);
    const double M = M0 + n * (t - t0);
    const double E = solve_kepler(M, k0.e);
    const double nu = ecsrt_v_ist(E, k0.e);
    const Kepler k{k0.a, k0.e, k0.inc, k0.raan, k0.argp, nu};
    return convert(k, mu);
}