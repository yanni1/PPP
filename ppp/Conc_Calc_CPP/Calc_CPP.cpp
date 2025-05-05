#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <vector>
#include <omp.h>
#include <iostream>
#include "Params.h"
#include <nanobind/stl/vector.h>
#include "Conc_Calc_CPP.cpp"

vector<vector<vector<float> > > calc(int nt_given) {
    //init C matrix (all values set to 0)
    params p = params(nt_given);
    vector<vector<vector<float> > > C(3, vector<vector<float> > (p.ny, vector<float>(p.nx, 0.0))); // C will contain C_CO C_CO2 and C_O2 => C[species,[y,[x,conc_value]]]
    //C[0] = C_CO
    //C[1] = C_CO2
    //C[2] = C_O2
    // CO concetration at y = 0 => C0_reservoir
    fill(C[0][0].begin(), C[0][0].end(), p.CO_reservoir);
    for (int n = 0; n < p.nt; n++){
        C = CC_CPP(nt_given, C);
    };
    
    return C;
}

//make nb module
NB_MODULE(calc, m) {
    m.def("calc", &calc, "calculates the conc gradient");
};