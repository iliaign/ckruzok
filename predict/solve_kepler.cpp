// solve_kepler.cpp
#include "solve_kepler.hpp"

#include <cmath>
#include <limits>

static constexpr int    kMaxIter = 50;
static constexpr double kTol     = std::numeric_limits<double>::epsilon() * 100.0;  // ~2.2e-14

double solve_kepler(const double M, const double e) {
    const double M_wrapped = std::remainder(M, 2.0 * M_PI);
    const double m = std::fabs(M_wrapped);
    double E = (e < 0.8) ? m : M_PI;

    for (int i = 0; i < kMaxIter; ++i) {
        const double dE = (E - e * std::sin(E) - m) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::fabs(dE) < kTol) break;
    }

    return std::copysign(E, M_wrapped);
}