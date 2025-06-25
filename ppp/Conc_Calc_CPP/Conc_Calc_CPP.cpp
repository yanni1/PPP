/*
 *  C++ source file for module ppp.Conc_Calc_CPP
 */


// See http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html for examples on how to use pybind11.
// The example below is modified after http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#More-on-working-with-numpy-arrays

#ifndef PROFILING
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <nanobind/stl/vector.h>
namespace nb = nanobind;
#endif

#include <vector>
#include <omp.h>
#include <iostream>
#include "Params.h"
#include "Reac_Rate_CPP.h"

using namespace std;

void CC_CPP(const params& p, vector<float>& C, vector<float>& Cn) { //make the passed vectors pointers (buffers)

    #pragma omp parallel for collapse(2)
    for (int j = 1; j < p.ny-1; j++){
        for (int i = 1; i < p.nx-1; i++){
            // circular O2 vent
            int x_dist2 = (i-p.a) * (i-p.a);
            int y_dist2 = (j-p.b) * (j-p.b);
            int distSq = x_dist2 + y_dist2;
            if (distSq <= p.r2){
                int o2_idx = p.idx(2, j, i);
                C[o2_idx] = p.O2_reservoir;
            };

            //diffusion
            float D2x_CO = float((C[p.idx(0,j,i+1)] + C[p.idx(0,j,i-1)] - 2*C[p.idx(0,j,i)]) / (p.dx * p.dx));
            float D2y_CO = float((C[p.idx(0,j+1,i)] + C[p.idx(0,j-1,i)] - 2*C[p.idx(0,j,i)]) / (p.dy * p.dy));
            float Diff_CO = p.D * (D2x_CO + D2y_CO);

            float D2x_CO2 = float((C[p.idx(1,j,i+1)] + C[p.idx(1,j,i-1)] - 2*C[p.idx(1,j,i)]) / (p.dx * p.dx));
            float D2y_CO2 = float((C[p.idx(1,j+1,i)] + C[p.idx(1,j-1,i)] - 2*C[p.idx(1,j,i)]) / (p.dy * p.dy));
            float Diff_CO2 = p.D * (D2x_CO2 + D2y_CO2);

            float D2x_O2 = float((C[p.idx(2,j,i+1)] + C[p.idx(2,j,i-1)] - 2*C[p.idx(2,j,i)]) / (p.dx * p.dx));
            float D2y_O2 = float((C[p.idx(2,j+1,i)] + C[p.idx(2,j-1,i)] - 2*C[p.idx(2,j,i)]) / (p.dy * p.dy));
            float Diff_O2 = p.D * (D2x_O2 + D2y_O2);

            //advection
            float dy_CO = float((C[p.idx(0,j,i)] - C[p.idx(0,j-1,i)]) / p.dy);
            float advec_CO = p.v[i] * dy_CO;

            float dy_CO2 = float((C[p.idx(1,j,i)] - C[p.idx(1,j-1,i)]) / p.dy);
            float advec_CO2 = p.v[i] * dy_CO2;

            float dy_O2 = float((C[p.idx(2,j,i)] - C[p.idx(2,j-1,i)])  / p.dy);
            float advec_O2 = p.v[i] * dy_O2;

            //reaction
            float k = reaction_rate_cpp(C, j, i, p);
            float reac_CO = k * C[p.idx(0,j,i)];
            float reac_CO2 = -1 * float(reac_CO); //right side of reaction
            float reac_O2 = float(reac_CO);

            //carbon capture ellipse
            float eps_CO = 0;
            float eps_CO2 = 0;
            if (((((i-p.x0_cc)*(i-p.x0_cc))/(p.semiMaj*p.semiMaj))+(((j-p.y0_cc)*(j-p.y0_cc))/(p.semiMin*p.semiMin))) <= 1.00f){
                float sum2_C = ((C[p.idx(0, j, i)] - p.mu_CO)*(C[p.idx(0, j, i)] - p.mu_CO)) + ((C[p.idx(1, j, i)] - p.mu_CO2)*(C[p.idx(1, j, i)] - p.mu_CO2));
                float rat_CO = (C[p.idx(0, j, i)] / sum2_C);
                float rat_CO2 = (C[p.idx(1, j, i)] / sum2_C);
                eps_CO = (p.alpha /(p.sigma_CO * sqrtf(2.0f * M_PI))) * expf(-0.5f*(rat_CO/(p.sigma_CO*p.sigma_CO)));
                eps_CO2 = (p.alpha /(p.sigma_CO2 * sqrtf(2.0f * M_PI))) * expf(-0.5f*(rat_CO2/(p.sigma_CO2*p.sigma_CO2)));
            };

            //update conc
            Cn[p.idx(0,j,i)] = C[p.idx(0,j,i)] + (p.dt * (Diff_CO - advec_CO - reac_CO)) - (eps_CO * C[p.idx(0,j,i)]);
            Cn[p.idx(1,j,i)] = C[p.idx(1,j,i)] + (p.dt * (Diff_CO2 - advec_CO2 - reac_CO2)) - (eps_CO2 * C[p.idx(1,j,i)]);
            Cn[p.idx(2,j,i)] = C[p.idx(2,j,i)] + (p.dt * (Diff_O2 - advec_O2 - reac_O2));

            
        };
    };
    //edges
    for (int j = 0; j < p.ny; j++) {
        Cn[p.idx(0,j,0)] = Cn[p.idx(0,j,1)]; //setting vertical edge (left) conc to the same as column next to it
        Cn[p.idx(1,j,0)] = Cn[p.idx(1,j,1)];
        Cn[p.idx(2,j,0)] = Cn[p.idx(2,j,1)];
        Cn[p.idx(0,j,p.nx - 1)] = Cn[p.idx(0,j,p.nx - 2)]; //setting vertical edge (right) conc to the same as column next to it
        Cn[p.idx(1,j,p.nx - 1)] = Cn[p.idx(1,j,p.nx - 2)]; 
        Cn[p.idx(2,j,p.nx - 1)] = Cn[p.idx(2,j,p.nx - 2)]; 
    }
    //re set CO boundary condition
    fill(C.begin() + p.idx(0, 0, 0), C.begin() + p.idx(0, 0, 0) + p.nx, p.CO_reservoir); //fill O's from 0 to 0+ nx (first row of species 0 (CO))
    return;
}

#ifndef PROFILING
//make nb module
NB_MODULE(Conc_Calc_CPP, m) {
    m.def("CC_CPP", &CC_CPP, "updates conc each timestep");
};
#endif
