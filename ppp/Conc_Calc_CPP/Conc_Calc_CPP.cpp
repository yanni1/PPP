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

void CC_CPP(const params& p, vector<float>& C, vector<float>& Cn, vector<float>& eps_field, bool& update_O2) { //make the passed vectors pointers (buffers)
    ofstream log_file("O2vent_log.txt");
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < p.ny-1; j++){
        for (int i = 1; i < p.nx-1; i++){

            // circular O2 vent
            if (update_O2 == true){
                log_file << "update_O2 was true and C[o2_idx] is now: ";
                int x_dist2 = (i-p.a) * (i-p.a);
                int y_dist2 = (j-p.b) * (j-p.b);
                int distSq = x_dist2 + y_dist2;
                if (distSq <= p.r2){
                    int o2_idx = p.idx(2, j, i);
                    C[o2_idx] += p.rho;
                    log_file << C[o2_idx]<<'\n';

                };
            }
            

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
            float eps_CO_raw = 0.0f;
            float eps_CO2_raw = 0.0f;
            float eps_CO = 0.0f;
            float eps_CO2 = 0.0f;
            int Ellips_Cond = (((i-p.x0_cc)*(i-p.x0_cc))/(p.semiMaj*p.semiMaj))+(((j-p.y0_cc)*(j-p.y0_cc))/(p.semiMin*p.semiMin));
            if (Ellips_Cond < 1){
                float sum_C = C[p.idx(0, j, i)] + C[p.idx(1, j, i)];
                if (sum_C > 1e-8f) { //add div by 0 protection => fix weird behaviour
                    float C_norm_CO = C[p.idx(0, j, i)]/sum_C;
                    float C_norm_CO2 = C[p.idx(1, j, i)]/sum_C;
                    float exp_CO = ((C_norm_CO - p.mu_CO) / p.sigma_CO) * ((C_norm_CO - p.mu_CO) / p.sigma_CO);
                    float exp_CO2 = ((C_norm_CO2 - p.mu_CO2) / p.sigma_CO2) * ((C_norm_CO2 - p.mu_CO2) / p.sigma_CO2);

                    //original square eps_field ellipse (should not be possible to be square)
                    eps_CO_raw = (p.alpha /(p.sigma_CO * sqrtf(2.0f * M_PI))) * expf(-0.5f*exp_CO);
                    eps_CO2_raw = (p.alpha /(p.sigma_CO2 * sqrtf(2.0f * M_PI))) * expf(-0.5f*exp_CO2);

                    //elliptical shape (converting my square pain to a nice ellipse, might be cheating)
                    float dx_norm = float(i - p.x0_cc) / float(p.semiMaj);
                    float dy_norm = float(j - p.y0_cc) / float(p.semiMin);
                    float norm_r2 = dx_norm * dx_norm + dy_norm * dy_norm;

                    float ellipse_falloff = 1.0f - norm_r2; // goes from 1 at center to 0 at edge
                    ellipse_falloff = fmaxf(0.0f, ellipse_falloff); // clip to [0,1]
                    
                    //ellipse conversion
                    eps_CO = ellipse_falloff * eps_CO_raw;
                    eps_CO2 = ellipse_falloff * eps_CO2_raw;
                    eps_field[j*p.nx + i] = eps_CO;

                } else { //if div by 0 safeguard is not met
                    eps_CO = 0;
                    eps_CO2 = 0;
                    eps_field[j*p.nx + i] = 0; //zeroeing eps to prevent "stale" values => data left in mem from previous steps
                };

            } else { //if ellipse cond is not met
                eps_CO = 0;
                eps_CO2 = 0;
                eps_field[j*p.nx + i] = 0; 
            };

            //update conc
            Cn[p.idx(0,j,i)] = C[p.idx(0,j,i)] + (p.dt * (Diff_CO - advec_CO - reac_CO)) - (p.dt * eps_CO * C[p.idx(0,j,i)])*100;
            Cn[p.idx(1,j,i)] = C[p.idx(1,j,i)] + (p.dt * (Diff_CO2 - advec_CO2 - reac_CO2)) - (p.dt * eps_CO2 * C[p.idx(1,j,i)])*100;
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
