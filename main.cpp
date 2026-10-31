#include <cmath>
#include <iostream>

#include "dec_to_kep.hpp"
#include "kep_to_dec.hpp"

int main() {
    using namespace kep_dec;

    const double mu = 398600.4418;  // Земля

    // Список тестов: {название, a, e, i, raan, argp, nu} (углы в градусах)
    struct Test { const char* name; double a, e, i, raan, argp, nu; };
    const Test tests[] = {
        {"Эллипс",        7000, 0.1, 45, 30, 60, 90},
        {"Круговая",      7000, 0.0, 45, 30,  0, 90},
        {"Экваториал",    7000, 0.2,  0,  0, 60, 90},
        {"Круг+экватор",  7000, 0.0,  0,  0,  0, 90},
    };

    const double d2r = M_PI / 180.0;

    for (const auto& t : tests) {
        const KeplerIn kep{t.a, t.e, t.i * d2r, t.raan * d2r, t.argp * d2r, t.nu * d2r};

        const auto [r, v]    = convert(kep, mu);    // Kep -> Dec
        const KeplerOut out  = convert(r, v, mu);   // Dec -> Kep
        const KeplerIn kep2{out.a, out.e, out.inc, out.raan, out.argp, out.nu};
        const auto [r2, v2]  = convert(kep2, mu);   // Kep -> Dec (проверка)

        std::cout << t.name
                  << ": ошибка r = " << (r - r2).norm()
                  << ", ошибка v = " << (v - v2).norm() << "\n";
    }
    return 0;
}
