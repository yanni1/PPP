#ifndef PROFILING
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <nanobind/stl/vector.h>
#endif

#include <vector>
//#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"

vector<float> calc(int nt_given) {
    //params 1 keer init en dan altijd als const doorgeven
    params p = params(nt_given);
    cout << "nt set to: " << nt_given << '\n' <<     endl;
   
    //flat concentration vectors [nt, ns, ny, nx]
    vector<float> Ct(nt_given * p.ns * p.ny * p.nx, 0.0f); //flattened vector values set to 0f (index = (s * ny + y) * nx + x)
    vector<float> C(p.ns * p.ny * p.nx, 0.0f); //flattened input buffer
    vector<float> Cn(p.ns * p.ny * p.nx, 0.0f); //flattened output buffer
    //C[0] = C_CO
    //C[1] = C_CO2
    //C[2] = C_O2

    // CO concetration at y = 0 => C0_reservoir
    //flattened initial condition
    fill(C.begin() + p.idx(0, 0, 0), C.begin() + p.idx(0, 0, 0) + p.nx, p.CO_reservoir);

    for (int n = 0; n < nt_given; n++) {
        cout << "timestep:" << n << '\n' << endl;
        CC_CPP(p, C, Cn);  // evolve concentration directly in preallocated vectors, only passed as pointers => using 2 buffers essentially
        //Store current timestep C into flat Ct
        for (int s = 0; s < p.ns; ++s) {
            for (int y = 0; y < p.ny; ++y) {
                for (int x = 0; x < p.nx; ++x) {
                    int ct_idx = ((n * p.ns + s) * p.ny + y) * p.nx + x;
                    int c_idx = (s * p.ny + y) * p.nx + x;
                    Ct[ct_idx] = C[c_idx];
                }
            }
        }
        swap(C, Cn); // reuse buffers for next step without copying => C@n -> C@n+1 (input) and Cn now holds old state (C@n) => does not matter, gets replaced by computed data anyway
    };                
    return Ct;
}

#ifndef PROFILING
//make nb module
NB_MODULE(Calc_CPP, m) {
    m.def("Calc_CPP", &calc, "calculates the conc gradient");
};
#endif