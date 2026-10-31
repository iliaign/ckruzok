#include "predict.hpp"

std::tuple<Eigen::Vector3d, Eigen::Vector3d> predict(const Eigen::Vector3d& r0,
                                                     const Eigen::Vector3d& v0,
                                                     const double t0,
                                                     const double t,
                                                     const double mu) {
    const Kepler k0 = convert(r0, v0, mu);               // 1) (r0, v0) -> элементы
    const double E0 = ist_v_ecrt(k0.nu, k0.e);           // 2) nu0 -> E0
    const double M0 = ean_mean(E0, k0.e);                // 3) E0 -> M0
    const double n  = mean_motion(k0.a, mu);             // 4) среднее движение
    const double M  = M0 + n * (t - t0);                 //    M = M0 + n (t - t0)
    const double E  = solve_kepler(M, k0.e);             // 5) M -> E
    const double nu = ecsrt_v_ist(E, k0.e);              // 6) E -> nu
    const Kepler k{k0.a, k0.e, k0.inc, k0.raan, k0.argp, nu};
    return convert(k, mu);                               // 7) элементы -> (r, v)
}