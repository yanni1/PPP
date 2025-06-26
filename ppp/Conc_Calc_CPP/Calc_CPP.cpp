#ifndef PROFILING
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h> 
#endif

#include <vector>
#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"

tuple<vector<float>, vector<float>>calc(int nt_given) {
    //params 1 keer init en dan altijd als const doorgeven
    params p = params(nt_given);
    cout << "nt set to: " << nt_given << '\n' <<     endl;
   
    //flat concentration vectors [nt, ns, ny, nx]
    vector<float> Ct(nt_given * p.ns * p.ny * p.nx, 0.0f); //flattened vector values set to 0f (index = (s * ny + y) * nx + x)
    vector<float> C(p.ns * p.ny * p.nx, 0.0f); //flattened input buffer
    vector<float> Cn(p.ns * p.ny * p.nx, 0.0f); //flattened output buffer

    //dubug
    vector<float> eps_field(p.ny * p.nx, 0.0f);
    vector<float> eps_series(nt_given * p.ny * p.nx, 0.0f);
    int spatial_size = p.ny * p.nx;
    //C[0] = C_CO
    //C[1] = C_CO2
    //C[2] = C_O2

    // CO concetration at y = 0 => C0_reservoir
    //flattened initial condition
    fill(C.begin() + p.idx(0, 0, 0), C.begin() + p.idx(0, 0, 0) + p.nx, p.CO_reservoir);

    for (int n = 0; n < nt_given; n++) {
        if (n % 100){cout << "timestep:" << n << '\n' << endl;};
        CC_CPP(p, C, Cn, eps_field);  // evolve concentration directly in preallocated vectors, only passed as pointers => using 2 buffers essentially
        //Store current timestep C into flat Ct
        int t_offset = n * spatial_size;
        #pragma omp parallel for simd
        for (int idx = 0; idx < spatial_size; idx++) {
            eps_series[t_offset + idx] = eps_field[idx];
        };
        #pragma omp parallel for simd
        for (int idx = 0; idx < p.ns * p.ny * p.nx; idx++) { //fully flattened data copying
            Ct[n * p.ns * p.ny * p.nx + idx] = C[idx];
        };
        swap(C, Cn); // reuse buffers for next step without copying => C@n -> C@n+1 (input) and Cn now holds old state (C@n) => does not matter, gets replaced by computed data anyway
    };                
    return {Ct, eps_series};
}

#ifndef PROFILING
//make nb module
//NB_MODULE(Calc_CPP, m) {
    //m.def("Calc_CPP", &calc, "calculates the conc gradient");
//};
NB_MODULE(Calc_CPP, m) {
    m.def("Calc_CPP", &calc, "calculates the conc gradient and eps_map");
};
#endif