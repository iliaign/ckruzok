#include <iostream>
#include "dec_to_kep.hpp"
#include "kep_to_dec.hpp"

int main() {
    double mu = 398600.4418; // Земля
    DecToKep d2k(mu);
    KepToDec k2d(mu);

    // Список тестов: {название, a, e, i, raan, argp, nu} (углы в градусах)
    struct Test { const char* name; double a, e, i, raan, argp, nu; };
    Test tests[] = {
        {"Эллипс",        7000, 0.1, 45, 30, 60, 90},
        {"Круговая",      7000, 0.0, 45, 30,  0, 90},
        {"Экваториал",    7000, 0.2,  0,  0, 60, 90},
        {"Круг+экватор",  7000, 0.0,  0,  0,  0, 90},
    };

    for (auto& t : tests) {
        Eigen::Matrix<double, 6, 1> kep;
        kep << t.a, t.e, t.i*M_PI/180, t.raan*M_PI/180, t.argp*M_PI/180, t.nu*M_PI/180;

        Eigen::Vector3d r, v, r2, v2;
        k2d.convert(kep, r, v);       // Kep -> Dec
        auto kep2 = d2k.convert(r, v); // Dec -> Kep
        k2d.convert(kep2, r2, v2);     // Kep -> Dec (проверка)

        std::cout << t.name << ": ошибка r = " << (r - r2).norm() << "\n";
    }
    return 0;
}