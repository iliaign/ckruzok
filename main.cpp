#include <cmath>
#include <iostream>

#include "kep_dec.hpp"
#include "predict.hpp"

int main() {

    const double mu = 398600.4418;  // Земля
    const double d2r = M_PI / 180.0;

    // Список тестов: {название, a, e, i, raan, argp, nu, число периодов}
    // (углы в градусах)
    struct Test { const char* name; double a, e, i, raan, argp, nu, periods; };
    const Test tests[] = {
        {"Эллипс",        7000, 0.1, 45, 30, 60, 90, 1.0},
        {"Эллипс 3 витка",7000, 0.1, 45, 30, 60, 90, 3.0},
        {"Большое e",     9000, 0.9, 45, 30, 60, 90, 1.0},
        {"Круговая",      7000, 0.0, 45, 30,  0, 90, 1.0},
        {"Экваториал",    7000, 0.2,  0,  0, 60, 90, 1.0},
    };

    for (const auto& t : tests) {
        const Kepler kep{t.a, t.e, t.i * d2r, t.raan * d2r, t.argp * d2r, t.nu * d2r};

        const double E_a  = ist_v_ecrt(kep.nu, kep.e);
        const double nu_a = ecsrt_v_ist(E_a, kep.e);
        const double M_a  = ean_mean(E_a, kep.e);
        const double E_b  = solve_kepler(M_a, kep.e);

        std::cout << t.name
                  << ": ошибка nu = " << std::fabs(std::remainder(nu_a - kep.nu, 2.0 * M_PI))
                  << ", ошибка E = "  << std::fabs(E_b - E_a);

        const auto [r0, v0] = convert(kep, mu);             
        const double dt = t.periods * 2.0 * M_PI / mean_motion(kep.a, mu);
        const auto [r1, v1] = predict(r0, v0, 0.0, dt, mu);  

        std::cout << ", ошибка r = " << (r1 - r0).norm()
                  << ", ошибка v = " << (v1 - v0).norm() << "\n";
    }
    return 0;
}