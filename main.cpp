// main.cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

#include "kep_dec/dec_to_kep.hpp"
#include "kep_dec/kep_to_dec.hpp"

const double MU = 398600.4418;
const double PI = 3.14159265358979323846;
const double DEG2RAD = PI / 180.0;
const double RAD2DEG = 180.0 / PI;

void header(const std::string& name) {
    std::cout << "\n============================================================\n";
    std::cout << "  " << name << "\n";
    std::cout << "============================================================\n";
}

void test_from_state(const std::string& name,
                     const Eigen::Vector3d& r, const Eigen::Vector3d& v)
{
    header(name);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  r = [" << r.transpose() << "] км\n";
    std::cout << "  v = [" << v.transpose() << "] км/с\n\n";

    Eigen::Matrix<double, 6, 1> kep = dec_to_kep(r, v, MU);

    std::cout << "  a    = " << kep(0) << " км\n";
    std::cout << "  e    = " << kep(1) << "\n";
    std::cout << "  i    = " << kep(2) * RAD2DEG << " град\n";
    std::cout << "  raan = " << kep(3) * RAD2DEG << " град\n";
    std::cout << "  argp = " << kep(4) * RAD2DEG << " град\n";
    std::cout << "  nu   = " << kep(5) * RAD2DEG << " град\n\n";

    Eigen::Vector3d r2, v2;
    kep_to_dec(kep, MU, r2, v2);

    double dr = (r2 - r).norm();
    double dv = (v2 - v).norm();

    std::cout << "  r2 = [" << r2.transpose() << "] км\n";
    std::cout << "  v2 = [" << v2.transpose() << "] км/с\n";
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  |r2 - r| = " << dr << " км\n";
    std::cout << "  |v2 - v| = " << dv << " км/с\n";

    bool ok = (dr < 1e-6) && (dv < 1e-9);
    std::cout << (ok ? "  [OK]\n" : "  [FAIL]\n");
}

int main() {
    {
        double r = 6778.0;
        double vc = std::sqrt(MU / r);
        Eigen::Vector3d rr(r, 0, 0);
        Eigen::Vector3d vv(0, vc * std::cos(51.6 * DEG2RAD),
                              vc * std::sin(51.6 * DEG2RAD));
        test_from_state("Круговая #1: LEO, i = 51.6", rr, vv);
    }

    {
        double r = 7000.0;
        double vc = std::sqrt(MU / r);
        Eigen::Vector3d rr(r, 0, 0);
        Eigen::Vector3d vv(0, 0, vc);
        test_from_state("Круговая #2: полярная, i = 90", rr, vv);
    }

    {
        double rp = 6600.0, ra = 42000.0;
        double a  = (rp + ra) / 2.0;
        double vp = std::sqrt(MU * (2.0 / rp - 1.0 / a));
        double inc = 28.5 * DEG2RAD;
        Eigen::Vector3d rr(rp, 0, 0);
        Eigen::Vector3d vv(0, vp * std::cos(inc), vp * std::sin(inc));
        test_from_state("Эллипс #1: GTO, i = 28.5", rr, vv);
    }

    {
        double rp = 6900.0, ra = 12800.0;
        double a  = (rp + ra) / 2.0;
        double vp = std::sqrt(MU * (2.0 / rp - 1.0 / a));
        double inc = 63.4 * DEG2RAD;
        Eigen::Vector3d rr(rp, 0, 0);
        Eigen::Vector3d vv(0, vp * std::cos(inc), vp * std::sin(inc));
        test_from_state("Эллипс #2: i = 63.4", rr, vv);
    }

    {
        double rp = 7000.0;
        double vp = std::sqrt(2.0 * MU / rp);
        Eigen::Vector3d rr(rp, 0, 0);
        Eigen::Vector3d vv(0, vp, 0);
        test_from_state("Парабола #1: rp = 7000, i = 0", rr, vv);
    }

    {
        double rp = 8000.0;
        double vp = std::sqrt(2.0 * MU / rp);
        double inc = 45.0 * DEG2RAD;
        Eigen::Vector3d rr(rp, 0, 0);
        Eigen::Vector3d vv(0, vp * std::cos(inc), vp * std::sin(inc));
        test_from_state("Парабола #2: rp = 8000, i = 45", rr, vv);
    }

    return 0;
}