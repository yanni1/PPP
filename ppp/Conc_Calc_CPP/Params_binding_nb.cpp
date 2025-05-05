#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include "Params.h"

namespace nb = nanobind;

NB_MODULE(params_module, m) {
    nb::class_<params>(m, "params")
        .def(nb::init<int>())
        .def_readwrite("dx", &params::dx)
        .def_readwrite("dy", &params::dy)
        .def_readwrite("dt", &params::dt)
        .def_readwrite("Lx", &params::Lx)
        .def_readwrite("Ly", &params::Ly)
        .def_readwrite("D", &params::D)
        .def_readwrite("k0", &params::k0)
        .def_readwrite("vmax", &params::vmax)
        .def_readwrite("CO_reservoir", &params::CO_reservoir)
        .def_readwrite("O2_reservoir", &params::O2_reservoir)
        .def_readwrite("nx", &params::nx)
        .def_readwrite("ny", &params::ny)
        .def_readwrite("nt", &params::nt)
        .def_readwrite("a", &params::a)
        .def_readwrite("b", &params::b)
        .def_readwrite("r", &params::r)
        .def_readwrite("r2", &params::r2)
        .def_readwrite("x", &params::x)
        .def_readwrite("y", &params::y)
        .def_readwrite("v", &params::v);
}
