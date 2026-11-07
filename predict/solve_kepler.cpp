// solve_kepler.cpp
#include "solve_kepler.hpp"

#include <cmath>



static constexpr int    kMaxIter = 50;
//тут подключить стд лимится
static constexpr double kTol     = 1e-14;  // остановка по |dE|



double solve_kepler(const double M, const double e) {
    // M -> [-pi, pi]; решаем для |M|, знак возвращаем в конце (уравнение нечётное)
    const double M_wrapped = std::remainder(M, 2.0 * M_PI);
    const double m = std::fabs(M_wrapped);

    // начальное приближение: для больших e стартуем с pi,
    // тогда ньютон сходится монотонно (f выпукла на [0, pi])
    double E = (e < 0.8) ? m : M_PI;

    for (int i = 0; i < kMaxIter; ++i) {
        const double dE = (E - e * std::sin(E) - m) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::fabs(dE) < kTol) break;
    }

    return std::copysign(E, M_wrapped);
}
