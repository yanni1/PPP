#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include "Params.h"

namespace nb = nanobind;

NB_MODULE(params_module, m) {
    nb::class_<params>(m, "params")
        .def(nb::init<int>())
        .def_prop_rw("dx", [](params &self) { return self.dx; }, [](params &self, float value) { self.dx = value; })
        .def_prop_rw("dy", [](params &self) { return self.dy; }, [](params &self, float value) { self.dy = value; })
        .def_prop_rw("dt", [](params &self) { return self.dt; }, [](params &self, float value) { self.dt = value; })
        .def_prop_rw("Lx", [](params &self) { return self.Lx; }, [](params &self, float value) { self.Lx = value; })
        .def_prop_rw("Ly", [](params &self) { return self.Ly; }, [](params &self, float value) { self.Ly = value; })
        .def_prop_rw("D", [](params &self) { return self.D; }, [](params &self, float value) { self.D = value; })
        .def_prop_rw("k0", [](params &self) { return self.k0; }, [](params &self, float value) { self.k0 = value; })
        .def_prop_rw("vmax", [](params &self) { return self.vmax; }, [](params &self, float value) { self.vmax = value; })
        .def_prop_rw("CO_reservoir", [](params &self) { return self.CO_reservoir; }, [](params &self, float value) { self.CO_reservoir = value; })
        .def_prop_rw("O2_reservoir", [](params &self) { return self.O2_reservoir; }, [](params &self, float value) { self.O2_reservoir = value; })
        .def_prop_rw("nx", [](params &self) { return self.nx; }, [](params &self, int value) { self.nx = value; })
        .def_prop_rw("ny", [](params &self) { return self.ny; }, [](params &self, int value) { self.ny = value; })
        .def_prop_rw("nt", [](params &self) { return self.nt; }, [](params &self, int value) { self.nt = value; })
        .def_prop_rw("a", [](params &self) { return self.a; }, [](params &self, float value) { self.a = value; })
        .def_prop_rw("b", [](params &self) { return self.b; }, [](params &self, float value) { self.b = value; })
        .def_prop_rw("r", [](params &self) { return self.r; }, [](params &self, float value) { self.r = value; })
        .def_prop_rw("r2", [](params &self) { return self.r2; }, [](params &self, float value) { self.r2 = value; })
        .def_prop_rw("x", [](params &self) { return self.x; }, [](params &self, std::vector<float> value) { self.x = value; })
        .def_prop_rw("y", [](params &self) { return self.y; }, [](params &self, std::vector<float> value) { self.y = value; })
        .def_prop_rw("v", [](params &self) { return self.v; }, [](params &self, std::vector<float> value) { self.v = value; });
}