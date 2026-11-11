#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include "dec_to_kep.hpp"
#include "ean_mean.hpp"
#include "ecsrt_v_ist.hpp"
#include "ist_v_ecrt.hpp"
#include "kep_to_dec.hpp"
#include "kepler.hpp"
#include "mean_motion.hpp"
#include "predict.hpp"
#include "solve_kepler.hpp"

namespace nb = nanobind;

NB_MODULE(pycosmos, m)
{
    m.doc() = "Orbit prediction library";

    nb::class_<Kepler>(m, "Kepler")
        .def(nb::init<>())
        .def(
            "__init__",
            [](Kepler *self, double a, double e, double inc, double raan, double argp, double nu) {
                new (self) Kepler{a, e, inc, raan, argp, nu};
            },
            nb::arg("a"), nb::arg("e"), nb::arg("inc"), nb::arg("raan"), nb::arg("argp"), nb::arg("nu"))
        .def_rw("a", &Kepler::a)
        .def_rw("e", &Kepler::e)
        .def_rw("inc", &Kepler::inc)
        .def_rw("raan", &Kepler::raan)
        .def_rw("argp", &Kepler::argp)
        .def_rw("nu", &Kepler::nu)
        .def("__repr__", [](const Kepler &k) {
            return "Kepler(a=" + std::to_string(k.a) + ", e=" + std::to_string(k.e) + ", inc=" + std::to_string(k.inc) +
                   ", raan=" + std::to_string(k.raan) + ", argp=" + std::to_string(k.argp) +
                   ", nu=" + std::to_string(k.nu) + ")";
        });

    // Две перегрузки convert -> два разных имени в Python
    m.def("kep_to_dec", nb::overload_cast<const Kepler &, double>(&convert), nb::arg("kep"), nb::arg("mu"),
          "Keplerian elements -> (r, v)");

    m.def("dec_to_kep", nb::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &, double>(&convert),
          nb::arg("r"), nb::arg("v"), nb::arg("mu"), "(r, v) -> Keplerian elements");

    m.def("predict", &predict, nb::arg("r0"), nb::arg("v0"), nb::arg("t0"), nb::arg("t"), nb::arg("mu"),
          "Propagate state from t0 to t; returns (r, v)");

    m.def("solve_kepler", &solve_kepler, nb::arg("M"), nb::arg("e"), nb::arg("max_iter") = 50);

    m.def("mean_motion", &mean_motion, nb::arg("a"), nb::arg("mu"));
    m.def("ist_v_ecrt", &ist_v_ecrt, nb::arg("nu"), nb::arg("e"));
    m.def("ecsrt_v_ist", &ecsrt_v_ist, nb::arg("E"), nb::arg("e"));
    m.def("ean_mean", &ean_mean, nb::arg("E"), nb::arg("e"));
}