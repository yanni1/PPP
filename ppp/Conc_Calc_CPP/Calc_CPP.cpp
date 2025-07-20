#ifndef PROFILING
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h> 
namespace nb = nanobind;
#endif

#include <vector>
#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Conc_Calc_CPP.h"

tuple<vector<float>, vector<float>>calc(params& p) {
    //params 1 keer init en dan altijd als const doorgeven

    cout << "nt set to: " << p.nt << '\n' << endl;
   
    //flat concentration vectors [nt, ns, ny, nx]
    vector<float> Ct(p.nt * p.ns * p.ny * p.nx, 0.0f); //flattened vector values set to 0f (index = (s * ny + y) * nx + x)
    vector<float> C(p.ns * p.ny * p.nx, 0.0f); //flattened input buffer
    vector<float> Cn(p.ns * p.ny * p.nx, 0.0f); //flattened output buffer
    
    //C[0] = C_CO
    //C[1] = C_CO2
    //C[2] = C_O2

    //absorption field visualization
    vector<float> eps_field(p.ny * p.nx, 0.0f);
    vector<float> eps_series(p.nt * p.ny * p.nx, 0.0f);
    int eps_size = p.ny * p.nx; //index eps_field

    // CO concetration at y = 0 => C0_reservoir
    //flattened initial condition
    fill(C.begin() + p.idx(0, 0, 0), C.begin() + p.idx(0, 0, 0) + p.nx, p.CO_reservoir);

    for (int n = 0; n < p.nt; n++) {
        if (n % 1000 == 0){cout << "timestep:" << n << '\n' << endl;};
        fill(eps_field.begin(), eps_field.end(), 0.0f); //reset eps_field
        //GA logic to start O2 vent and activate pump
        int update_O2 = 0;
        if ( n == 0){
            update_O2 = 1; //init O2 vent with reservoir
        } else if (n % p.tau == 0) {
            update_O2 = 2; //activate pump every tau time steps
        };

        CC_CPP(p, C, Cn, eps_field, update_O2);  // evolve concentration directly in preallocated vectors, only passed as pointers => using 2 buffers essentially
        
        //Store current timestep C into flat Ct and same for eps_field
        int t_offset = n * eps_size;
        #pragma omp parallel for simd
        for (int idx = 0; idx < eps_size; idx++) {
            eps_series[t_offset + idx] = eps_field[idx]; //eps buffer not swapped just reused (reset to 0)
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
NB_MODULE(Calc_CPP, m) {
    m.def("Calc_CPP", &calc, nb::rv_policy::move, "calculates the conc gradient and eps_map"); //using move to move the cpp vectors to python via nb instead of default copying
};
#endif