/*
 *  C++ source file for module ppp.Conc_Calc_CPP
 */

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
#include <fstream>

using namespace std;

void CC_CPP(const params& p, vector<float>& C, vector<float>& Cn, vector<float>& eps_field, int& update_O2) { //make the passed vectors pointers (buffers)
    //precomputed O2 checks
    const bool do_initial_O2 = (update_O2 == 1); //true if update_O2 = 1
    const bool do_pump_O2    = (update_O2 == 2); //true if update_O2 = 2
    // precomputed invariants
    const float inv_dx2 = 1.0f / (p.dx * p.dx);
    const float inv_dy2 = 1.0f / (p.dy * p.dy);
    const float inv_dy  = 1.0f / p.dy;
    //start main loop
    #pragma omp parallel for simd collapse(2)
    for (int j = 1; j < p.ny-1; j++){
        for (int i = 1; i < p.nx-1; i++){
            
            // O2 logic => compute distSq once per i and j update initial cond. and pump with precomputed bools
            if (do_initial_O2 || do_pump_O2) { //if is true als 1 van de twee true is
                int x_dist2 = (i - p.a) * (i - p.a);
                int y_dist2 = (j - p.b) * (j - p.b);
                int distSq  = x_dist2 + y_dist2;
                if (distSq <= p.r2){
                    int o2_idx = p.idx(2, j, i);
                    if (do_initial_O2){
                        C[o2_idx] = p.O2_reservoir;
                    } else if (do_pump_O2){
                        C[o2_idx] += p.rho;
                    };
                };
            };
            
            //diffusion
            float D2x_CO = float((C[p.idx(0,j,i+1)] + C[p.idx(0,j,i-1)] - 2*C[p.idx(0,j,i)]) * inv_dx2);
            float D2y_CO = float((C[p.idx(0,j+1,i)] + C[p.idx(0,j-1,i)] - 2*C[p.idx(0,j,i)]) * inv_dy2);
            float Diff_CO = p.D * (D2x_CO + D2y_CO);

            float D2x_CO2 = float((C[p.idx(1,j,i+1)] + C[p.idx(1,j,i-1)] - 2*C[p.idx(1,j,i)]) * inv_dx2);
            float D2y_CO2 = float((C[p.idx(1,j+1,i)] + C[p.idx(1,j-1,i)] - 2*C[p.idx(1,j,i)]) * inv_dy2);
            float Diff_CO2 = p.D * (D2x_CO2 + D2y_CO2);

            float D2x_O2 = float((C[p.idx(2,j,i+1)] + C[p.idx(2,j,i-1)] - 2*C[p.idx(2,j,i)]) * inv_dx2);
            float D2y_O2 = float((C[p.idx(2,j+1,i)] + C[p.idx(2,j-1,i)] - 2*C[p.idx(2,j,i)]) * inv_dy2);
            float Diff_O2 = p.D * (D2x_O2 + D2y_O2);

            //advection
            float dy_CO = float((C[p.idx(0,j,i)] - C[p.idx(0,j-1,i)]) * inv_dy);
            float advec_CO = p.v[i] * dy_CO;

            float dy_CO2 = float((C[p.idx(1,j,i)] - C[p.idx(1,j-1,i)]) * inv_dy);
            float advec_CO2 = p.v[i] * dy_CO2;

            float dy_O2 = float((C[p.idx(2,j,i)] - C[p.idx(2,j-1,i)])  * inv_dy);
            float advec_O2 = p.v[i] * dy_O2;

            //reaction
            float k = reaction_rate_cpp(C, j, i, p);
            float reac_CO = k * C[p.idx(0,j,i)];
            float reac_CO2 = -1 * float(reac_CO); //right side of reaction
            float reac_O2 = float(reac_CO);

            //carbon capture ellipse
            float eps_CO = 0.0f;
            float eps_CO2 = 0.0f;
            float Ellips_Cond = (static_cast<float>((i-p.x0_cc)*(i-p.x0_cc))/(p.semiMaj*p.semiMaj))+(static_cast<float>((j-p.y0_cc)*(j-p.y0_cc))/(p.semiMin*p.semiMin)); //cast floats to prevent int division from making a square
            if (Ellips_Cond <= 1.0f){
                float sum_C = C[p.idx(0, j, i)] + C[p.idx(1, j, i)];
                if (sum_C > 1e-8f) { //add div by 0 protection => fix weird behaviour
                    float C_norm_CO = C[p.idx(0, j, i)]/sum_C;
                    float C_norm_CO2 = C[p.idx(1, j, i)]/sum_C;
                    float exp_CO = ((C_norm_CO - p.mu_CO) / p.sigma_CO) * ((C_norm_CO - p.mu_CO) / p.sigma_CO);
                    float exp_CO2 = ((C_norm_CO2 - p.mu_CO2) / p.sigma_CO2) * ((C_norm_CO2 - p.mu_CO2) / p.sigma_CO2);

                    eps_CO = (p.alpha /(p.sigma_CO * sqrtf(2.0f * M_PI))) * expf(-0.5f*exp_CO);
                    eps_CO2 = (p.alpha /(p.sigma_CO2 * sqrtf(2.0f * M_PI))) * expf(-0.5f*exp_CO2);
                    eps_field[j * p.nx + i] = eps_CO2;   
                } 
            };

            //update conc
            Cn[p.idx(0,j,i)] = C[p.idx(0,j,i)] + (p.dt * (Diff_CO - advec_CO - reac_CO)) - (p.dt * eps_CO * C[p.idx(0,j,i)]);
            Cn[p.idx(1,j,i)] = C[p.idx(1,j,i)] + (p.dt * (Diff_CO2 - advec_CO2 - reac_CO2)) - (p.dt * eps_CO2 * C[p.idx(1,j,i)]);
            Cn[p.idx(2,j,i)] = C[p.idx(2,j,i)] + (p.dt * (Diff_O2 - advec_O2 - reac_O2)); 

            #ifdef PROFILING
            //safety checks //aborts if eps is not finite and positive
            assert(std::isfinite(eps_field[j * p.nx + i]));
            assert(eps_field[j * p.nx + i] >= 0.0f);
            #endif
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
    };

#ifndef PROFILING
//make nb module
NB_MODULE(Conc_Calc_CPP, m) {
    m.def("CC_CPP", &CC_CPP, "updates conc each timestep");
};
#endif
